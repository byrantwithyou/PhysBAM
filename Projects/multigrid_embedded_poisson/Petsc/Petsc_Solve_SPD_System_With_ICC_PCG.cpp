//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <cassert>
#include <cmath>

#include <algorithm>
#include <limits>
#include <ostream>
#include <vector>

#include <petsc.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscvec.h>

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <Jeffrey_Utilities/Algorithm/Fill.h>
#include <Jeffrey_Utilities/Algorithm/Find_All_If.h>
#include <Jeffrey_Utilities/Algorithm/For_Each.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/NOT_EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/GENERIC_SYSTEM_REFERENCE.h>
#include <Jeffrey_Utilities/Petsc/Add_Stencil_To_Matrix.h>
#include <Jeffrey_Utilities/Petsc/CALL_AND_CHKERRQ.h>
#include <Jeffrey_Utilities/Petsc/Print_KSP_Info.h>
#include <Jeffrey_Utilities/Petsc/SCOPED_DESTROY.h>
#include <Jeffrey_Utilities/POISSON_SOLVER_PARAMS.h>
#include <Jeffrey_Utilities/Stencils/SKIP_ZERO_VALUE_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL.h>
#include <Jeffrey_Utilities/Stencils/UNSTRUCTURED_STENCIL_PROXY.h>

#include <Jeffrey_Utilities/Petsc/Petsc_Solve_SPD_System_With_ICC_PCG.h>

namespace PhysBAM
{

namespace Petsc
{

namespace Detail_Solve_SPD_System_With_ICC_PCG
{

struct SET_PETSC_INDEX_OF_LINEAR_INDEX;

} // namespace Detail_Solve_SPD_System_With_ICC_PCG

template< class T >
PetscErrorCode
Solve_SPD_System_With_ICC_PCG(
    const unsigned int n_thread,
    const SOLVER_PARAMS params,
    const bool has_constant_vectors_in_null_space,
    const GENERIC_SYSTEM_REFERENCE<T> system,
    const ARRAY_VIEW<const T> rhs,
    ARRAY_VIEW<T> x,
    std::ostream& lout /*= PhysBAM::nout*/)
{
    typedef Detail_Solve_SPD_System_With_ICC_PCG::SET_PETSC_INDEX_OF_LINEAR_INDEX SET_PETSC_INDEX_OF_LINEAR_INDEX_;

    const int n_physbam_index = x.Size();
    assert(n_physbam_index == x.Size());
    assert(n_physbam_index == rhs.Size());

    BASIC_TIMER timer;
    MPI_Comm mpi_comm = PETSC_COMM_WORLD;

    lout << "Initializing PETSc index <-> PhysBAM index mappings...";
    lout.flush();
    timer.Restart();
    std::vector<int> physbam_index_of_petsc_index;
    Find_All_If_MT(
        n_thread,
        1, n_physbam_index,
        Make_Compose_Function(
            NOT_EQUAL1_FUNCTION<int>(0),
            PHYSBAM_BOUND_FAST_MEM_FN_TEMPLATE(
                system,
                &GENERIC_SYSTEM_REFERENCE<T>::Stencil_N_Nonzero
            )
        ),
        physbam_index_of_petsc_index
    );
    const int n_petsc = static_cast< int >(physbam_index_of_petsc_index.size());
    ARRAY<int> petsc_index_of_physbam_index(n_physbam_index, false); // uninit'ed
    Fill_MT(
        n_thread,
        petsc_index_of_physbam_index,
        std::numeric_limits<int>::max()
    );
    For_Each_MT(
        n_thread,
        0, n_petsc - 1,
        SET_PETSC_INDEX_OF_LINEAR_INDEX_(
            physbam_index_of_petsc_index,
            petsc_index_of_physbam_index
        )
    );
    lout << timer.Elapsed() << " s" << std::endl;
    lout << "  # of petsc dofs = " << n_petsc << std::endl;
    if(n_petsc == 0)
        return 0;

    lout << "Constructing Jacobi scaling coefficients...";
    lout.flush();
    timer.Restart();
    std::vector< double > jscalings(n_petsc); // init'ed to 0
    for(int i = 0; i != n_petsc; ++i) {
        const int physbam_index = physbam_index_of_petsc_index[i];
        assert(system.Stencil_N_Nonzero(physbam_index) != 0);
        const T diag = system.Diag(physbam_index);
        jscalings[i] = diag > 0 ? std::sqrt(static_cast< double >(diag)) : 1;
    }
    lout << timer.Elapsed() << " s" << std::endl;

    std::vector< PetscScalar > petsc_jns_vec_;
    Vec petsc_jns_vec;
    if(has_constant_vectors_in_null_space) {
        lout << "Constructing null-space vector of Jacobi-scaled matrix...";
        lout.flush();
        timer.Restart();
        petsc_jns_vec_.resize(n_petsc); // init'ed to 0
        // MatNullSpaceCreate requires the null-space vector to have unit 2-norm.
        const double norm = std::sqrt(std::inner_product(
            jscalings.begin(), jscalings.end(),
            jscalings.begin(),
            static_cast< double >(0)
        ));
        assert(norm > 0);
        for(int i = 0; i != n_petsc; ++i)
            petsc_jns_vec_[i] = static_cast< PetscScalar >(jscalings[i] / norm);
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( VecCreateSeqWithArray(
            mpi_comm,
            n_petsc,
            &petsc_jns_vec_.front(),
            &petsc_jns_vec
        ) );
    }
    PHYSBAM_PETSC_SCOPED_DESTROY_IF( Vec, petsc_jns_vec, has_constant_vectors_in_null_space );
    if(has_constant_vectors_in_null_space)
        lout << timer.Elapsed() << " s" << std::endl;

    lout << "Constructing matrix in CSR format and Jacobi scaling...";
    lout.flush();
    timer.Restart();
    // CSR = Compressed Sparse Row, aka Yale Sparse Matrix Format (according to Wikipedia)
    //       see http://en.wikipedia.org/wiki/Sparse_matrix#Yale_Sparse_Matrix_Format
    int matrix_nnz = 0;
    int max_stencil_nnz = 0;
    for(int i = 0; i != n_petsc; ++i) {
        const int physbam_index = physbam_index_of_petsc_index[i];
        const int stencil_nnz = system.Stencil_N_Nonzero(physbam_index);
        assert(stencil_nnz > 0);
        matrix_nnz += stencil_nnz;
        max_stencil_nnz = std::max(max_stencil_nnz, stencil_nnz);
    }
    std::vector< PetscInt > petsc_jmatrix_value_index_of_row(n_petsc + 1); // init'ed to 0
    std::vector< PetscInt > petsc_jmatrix_col_indices(matrix_nnz); // init'ed to 0
    std::vector< PetscScalar > petsc_jmatrix_values(matrix_nnz); // init'ed to 0
    {
        int value_index = 0;
        UNSTRUCTURED_STENCIL<int,T> system_stencil;
        system_stencil.values.Preallocate(max_stencil_nnz);
        UNSTRUCTURED_STENCIL_PROXY< UNSTRUCTURED_STENCIL<int,T> > system_stencil_proxy(system_stencil);
        for(int i = 0; i != n_petsc; ++i) {
            petsc_jmatrix_value_index_of_row[i] = value_index;
            const int physbam_index = physbam_index_of_petsc_index[i];
            system_stencil.values.Remove_All();
            system.Add_Stencil_To(physbam_index, system_stencil_proxy);
            const int next_value_index = Petsc::Add_Stencil_To_Matrix(
                petsc_index_of_physbam_index,
                Make_Skip_Zero_Value_Stencil_Proxy(system_stencil_proxy),
                value_index,
                petsc_jmatrix_col_indices,
                petsc_jmatrix_values
            );
            for(; value_index != next_value_index; ++value_index) {
                petsc_jmatrix_values[value_index] = static_cast< PetscScalar >(
                    petsc_jmatrix_values[value_index] /
                    (jscalings[i] * jscalings[petsc_jmatrix_col_indices[value_index]])
                );
            }
        }
        assert(value_index <= matrix_nnz);
        petsc_jmatrix_value_index_of_row.back() = value_index;
    }
    Mat petsc_jmatrix;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( MatCreateSeqAIJWithArrays(
        mpi_comm,
        n_petsc,
        n_petsc,
        &petsc_jmatrix_value_index_of_row.front(),
        &petsc_jmatrix_col_indices.front(),
        &petsc_jmatrix_values.front(),
        &petsc_jmatrix
    ) );
    PHYSBAM_PETSC_SCOPED_DESTROY( Mat, petsc_jmatrix );
    lout << timer.Elapsed() << " s" << std::endl;
    lout << "  # of nonzeros (estimated) = " << matrix_nnz << std::endl;
    lout << "  # of nonzeros (actual)    = " << petsc_jmatrix_value_index_of_row.back() << std::endl;

    // Print out how good this so-called null space vector is.
    if(has_constant_vectors_in_null_space) {
        std::vector< PetscScalar > petsc_v_(n_petsc); // init'ed to 0
        Vec petsc_v;
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( VecCreateSeqWithArray(mpi_comm, n_petsc, &petsc_v_.front(), &petsc_v) );
        PHYSBAM_PETSC_SCOPED_DESTROY( Vec, petsc_v );
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( MatMult(petsc_jmatrix, petsc_jns_vec, petsc_v) );
        PetscReal jnorm = static_cast< PetscReal >(0);
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( VecNorm(petsc_v, NORM_INFINITY, &jnorm) );
        lout << "  |(J^(-1)*A*J^(-1)) * (J*1)|_{infty} = " << jnorm << std::endl;
        PetscReal norm = 0;
        for(int i = 0; i != n_petsc; ++i)
            norm = std::max(norm, static_cast< PetscReal >(std::abs(petsc_v_[i] * jscalings[i])));
        lout << "  |J * (J^(-1)*A*J^(-1)) * (J*1)|_{infty} = " << norm << std::endl;
    }

    lout << "Constructing rhs of Jacobi-scaled system...";
    lout.flush();
    timer.Restart();
    std::vector< PetscScalar > petsc_jrhs_(n_petsc); // init'ed to 0
    for(int i = 0; i != n_petsc; ++i) {
        const int physbam_index = physbam_index_of_petsc_index[i];
        petsc_jrhs_[i] = static_cast< PetscScalar >(rhs(physbam_index) / jscalings[i]);
    }
    Vec petsc_jrhs;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( VecCreateSeqWithArray(mpi_comm, n_petsc, &petsc_jrhs_.front(), &petsc_jrhs) );
    PHYSBAM_PETSC_SCOPED_DESTROY( Vec, petsc_jrhs );
    lout << timer.Elapsed() << " s" << std::endl;

    // Construct the Petsc KSP object.
    KSP petsc_ksp;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPCreate(mpi_comm, &petsc_ksp) );
    PHYSBAM_PETSC_SCOPED_DESTROY( KSP, petsc_ksp );
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPSetType(petsc_ksp, KSPCG) );
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPSetOperators(petsc_ksp, petsc_jmatrix, petsc_jmatrix, SAME_NONZERO_PATTERN) );

    // Set the preconditioner as ICC (InComplete Cholesky).
    {
        PC petsc_pc;
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPGetPC(petsc_ksp, &petsc_pc) );
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( PCSetType(petsc_pc, params.precondition ? PCICC : PCNONE) );
    }

    // If has_constant_vectors_in_null_space, construct the null space (from the
    // null space vector) and let the CG solver know about it.
    MatNullSpace petsc_jnullspace;
    if(has_constant_vectors_in_null_space)
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( MatNullSpaceCreate(
            mpi_comm,
            PETSC_FALSE, // null space does *not* (generally) contain the constant vector
            1,           // dimension of null space
            &petsc_jns_vec,
            &petsc_jnullspace
        ) );
    PHYSBAM_PETSC_SCOPED_DESTROY_IF( MatNullSpace, petsc_jnullspace, has_constant_vectors_in_null_space );
    if(has_constant_vectors_in_null_space) {
        if(params.print_diagnostics) {
            PetscTruth is_null_space;
            PHYSBAM_PETSC_CALL_AND_CHKERRQ( MatNullSpaceTest(petsc_jnullspace, petsc_jmatrix, &is_null_space) );
            if(!is_null_space) {
                std::cerr << "ERROR: Petsc null space test failed." << std::endl;
                return 1;
            }
        }
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPSetNullSpace(petsc_ksp, petsc_jnullspace) );
    }

    // Initial guess will be rhs (which is generally nonzero) so that we
    // automatically satisfy any Dirichlet boundary conditions on the grid
    // boundary.
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPSetInitialGuessNonzero(petsc_ksp, PETSC_TRUE) );

    // Set additional solver options.
    if(params.print_residuals)
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPMonitorSet(petsc_ksp, &KSPMonitorTrueResidualNorm, PETSC_NULL, PETSC_NULL) );
    const unsigned int max_iterations =
        params.max_iterations == std::numeric_limits< unsigned int >::max() &&
        static_cast< unsigned int >(n_petsc) <= std::numeric_limits< unsigned int >::max() / 4 ?
        4 * static_cast< unsigned int >(n_petsc) :
        params.max_iterations;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPSetTolerances(
        petsc_ksp,
        static_cast< PetscReal >(params.relative_tolerance),
        static_cast< PetscReal >(params.absolute_tolerance),
        PETSC_DEFAULT, // divergence tolerance
        max_iterations
    ) );
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPSetComputeSingularValues(petsc_ksp, PETSC_TRUE) );
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPSetUp(petsc_ksp) );

    lout << "Executing KSPSolve..." << std::endl;
    timer.Restart();
    // petsc_jrhs will be overwritten with the solution.
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPSolve(petsc_ksp, petsc_jrhs, petsc_jrhs) );
    lout << "[Executing KSPSolve...] " << timer.Elapsed() << " s" << std::endl;

    if(params.print_diagnostics)
        PHYSBAM_PETSC_CALL_AND_CHKERRQ( KSPView(petsc_ksp, PETSC_VIEWER_STDOUT_SELF) );

    lout << "KSP info:" << std::endl;
    PHYSBAM_PETSC_CALL_AND_CHKERRQ( Print_KSP_Info(petsc_ksp, lout) );

    // Un-Jacobi scale the CG solution to get the actual approximate solution.
    for(int i = 0; i != n_petsc; ++i) {
        const int physbam_index = physbam_index_of_petsc_index[i];
        x(physbam_index) = static_cast<T>(petsc_jrhs_[i] / jscalings[i]);
    }

    return 0;
}

namespace Detail_Solve_SPD_System_With_ICC_PCG
{

struct SET_PETSC_INDEX_OF_LINEAR_INDEX
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        SET_PETSC_INDEX_OF_LINEAR_INDEX,
        (( std::vector<int> const &, physbam_index_of_petsc_index ))
        (( ARRAY<int>&, petsc_index_of_physbam_index ))
    )
public:
    typedef void result_type;
    void operator()(const int petsc_index) const
    {
        assert(0 <= petsc_index && static_cast< unsigned int >(petsc_index) < physbam_index_of_petsc_index.size());
        const int physbam_index = physbam_index_of_petsc_index[petsc_index];
        petsc_index_of_physbam_index(physbam_index) = petsc_index;
    }
};

} // namespace Detail_Solve_SPD_System_With_ICC_PCG

#define EXPLICIT_INSTANTIATION( T ) \
template \
PetscErrorCode \
Solve_SPD_System_With_ICC_PCG<T>( \
    const unsigned int n_thread, \
    const SOLVER_PARAMS params, \
    const bool has_constant_vectors_in_null_space, \
    const GENERIC_SYSTEM_REFERENCE<T> system, \
    const ARRAY_VIEW<const T> rhs, \
    ARRAY_VIEW<T> x, \
    std::ostream& lout);
EXPLICIT_INSTANTIATION( float )
EXPLICIT_INSTANTIATION( double )
#undef EXPLICIT_INSTANTIATION

} // namespace Petsc

} // namespace PhysBAM
