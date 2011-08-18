//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#include <cassert>
#include <cmath>

#include <algorithm>
#include <limits>
#include <ostream>

#include <Jeffrey_Utilities/Algorithm/Find_All_If.h>
#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/As_Sparse_Matrix_Flat_NXN.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/Functional/ARGUMENT_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/FILL_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/NULL_FUNCTION.h>
#include <Jeffrey_Utilities/GENERIC_SYSTEM_REFERENCE.h>
#include <Jeffrey_Utilities/Krylov/DOT_PRODUCT_INNER_PRODUCT.h>
#include <Jeffrey_Utilities/Krylov/JACOBI_PRECONDITIONER.h>
#include <Jeffrey_Utilities/Krylov/KRYLOV_SYSTEM_COMPOSER.h>
#include <Jeffrey_Utilities/Krylov/KRYLOV_VECTOR_WRAPPER_ARRAY_VIEW.h>
#include <Jeffrey_Utilities/Krylov/MAXIMUM_MAGNITUDE_CONVERGENCE_NORM.h>
#include <Jeffrey_Utilities/SOLVER_PARAMS.h>
#include <Jeffrey_Utilities/VISITOR_SEQUENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>

#include <Jeffrey_Utilities/Krylov/Solve_SPD_System_With_ICC_PCG.h>

namespace PhysBAM
{

template< class T >
void
Solve_SPD_System_With_ICC_PCG(
    const unsigned int n_thread,
    const SOLVER_PARAMS params,
    const bool has_constant_vectors_in_null_space,
    const GENERIC_SYSTEM_REFERENCE<T> system,
    const ARRAY_VIEW<T> rhs,
    ARRAY_VIEW<T> x,
    std::ostream& lout /*= PhysBAM::nout*/)
{
    assert(rhs.Size() == x.Size());

    BASIC_TIMER timer;

    float initial_residual_norm;
    {
        ARRAY<T> t(x.Size()); // init'ed to 0
        system.Apply(x, t);
        initial_residual_norm = static_cast< float >(ARRAYS_COMPUTATIONS::Maximum_Magnitude(t -= rhs));
    }
    const float tolerance = std::max(params.relative_tolerance * initial_residual_norm, params.absolute_tolerance);
    const int max_iterations =
        params.max_iterations == std::numeric_limits< unsigned int >::max()
     && x.Size() <= std::numeric_limits< int >::max() / 4
      ? 4 * x.Size()
      : static_cast< int >(std::min(
            params.max_iterations,
            static_cast< unsigned int >(std::numeric_limits< int >::max())
        ));

    if(params.precondition) {

        PCG_SPARSE<T> pcg;
        pcg.enforce_compatibility = has_constant_vectors_in_null_space;
        pcg.remove_null_space_solution_component = has_constant_vectors_in_null_space;
        pcg.show_residual = params.print_residuals;
        pcg.show_results = params.print_diagnostics;
        pcg.modified_incomplete_cholesky_coefficient = static_cast<T>(0.99);
        pcg.maximum_iterations = max_iterations;

        lout << "Constructing SPARSE_MATRIX_FLAT_NXN matrix...";
        lout.flush();
        timer.Restart();
        SPARSE_MATRIX_FLAT_NXN<T> flat_system;
        ARRAY<int> orig_index_of_flat_index;
        ARRAY<int> flat_index_of_orig_index(x.Size(), false); // uninit'ed
        As_Sparse_Matrix_Flat_NXN(
            n_thread,
            system,
            flat_system,
            orig_index_of_flat_index,
            flat_index_of_orig_index
        );
        lout << timer.Elapsed() << " s" << std::endl;
        const int n_flat_index = orig_index_of_flat_index.Size();
        lout << "  # of flat dofs = " << n_flat_index << std::endl;

        lout << "Constructing Jacobi scaling coefficients...";
        lout.flush();
        timer.Restart();
        ARRAY<double> jscalings(n_flat_index); // init'ed to 0
        for(int flat_index = 1; flat_index <= n_flat_index; ++flat_index) {
            const int orig_index = orig_index_of_flat_index(flat_index);
            assert(system.Stencil_N_Nonzero(orig_index) != 0);
            const T diag = system.Diag(orig_index);
            jscalings(flat_index) = diag > 0 ? std::sqrt(static_cast< double >(diag)) : 1;
        }
        lout << timer.Elapsed() << " s" << std::endl;

        lout << "Jacobi scaling SPARSE_MATRIX_FLAT_NXN matrix...";
        lout.flush();
        for(int flat_index = 1; flat_index <= n_flat_index; ++flat_index) {
            const int offset_index = flat_system.offsets(flat_index);
            const int next_offset_index = flat_system.offsets(flat_index + 1);
            for(int i = offset_index; i != next_offset_index; ++i) {
                const int other_flat_index = flat_system.A(i).j;
                flat_system.A(i).a = static_cast<T>(
                    flat_system.A(i).a / (jscalings(flat_index) * jscalings(other_flat_index))
                );
            }
        }
        lout << timer.Elapsed() << " s" << std::endl;

        lout << "Allocating storage for temporary vectors...";
        lout.flush();
        timer.Restart();
        VECTOR_ND<T> y(n_flat_index); // init'ed to 0
        VECTOR_ND<T> b(n_flat_index); // init'ed to 0
        VECTOR_ND<T> q(n_flat_index); // init'ed to 0
        VECTOR_ND<T> s(n_flat_index); // init'ed to 0
        VECTOR_ND<T> r(n_flat_index); // init'ed to 0
        VECTOR_ND<T> k(n_flat_index); // init'ed to 0
        VECTOR_ND<T> z(n_flat_index); // init'ed to 0
        for(int flat_index = 1; flat_index <= n_flat_index; ++flat_index) {
            const int orig_index = orig_index_of_flat_index(flat_index);
            y(flat_index) = static_cast<T>(x(orig_index) * jscalings(flat_index));
            b(flat_index) = static_cast<T>(rhs(orig_index) / jscalings(flat_index));
        }
        lout << timer.Elapsed() << " s" << std::endl;
        lout << "Executing PCG_SPARSE::Solve..." << std::endl;
        timer.Restart();
        pcg.Solve(
            flat_system,
            y, b, q, s, r, k, z,
            tolerance
        );
        lout << "[Executing PCG_SPARSE::Solve...] " << timer.Elapsed() << " s" << std::endl;

        for(int flat_index = 1; flat_index <= n_flat_index; ++flat_index) {
            const int orig_index = orig_index_of_flat_index(flat_index);
            x(orig_index) = static_cast<T>(y(flat_index) / jscalings(flat_index));
        }

    }
    else {

        // TODO: If non-trivial null space, use project_nullspace...?
        CONJUGATE_GRADIENT<T> cg;
        cg.print_diagnostics = params.print_diagnostics;
        cg.print_residuals = params.print_residuals;

        lout << "Allocating storage for temporary vectors...";
        lout.flush();
        timer.Restart();
        ARRAY<T> q(x.Size()); // init'ed to 0
        ARRAY<T> s(x.Size()); // init'ed to 0
        ARRAY<T> r(x.Size()); // init'ed to 0
        ARRAY<T> k(x.Size()); // init'ed to 0
        ARRAY<T> z(x.Size()); // init'ed to 0
        lout << timer.Elapsed() << " s" << std::endl;

        lout << "Executing CONJUGATE_GRADIENT::Solve..." << std::endl;
        timer.Restart();
        typedef KRYLOV_VECTOR_WRAPPER< T, ARRAY_VIEW<T> > KRYLOV_VECTOR_TYPE;
        KRYLOV_VECTOR_TYPE kx(x);
        KRYLOV_VECTOR_TYPE kb(rhs);
        KRYLOV_VECTOR_TYPE kq(As_Array_View(q));
        KRYLOV_VECTOR_TYPE ks(As_Array_View(s));
        KRYLOV_VECTOR_TYPE kr(As_Array_View(r));
        KRYLOV_VECTOR_TYPE kk(As_Array_View(k));
        KRYLOV_VECTOR_TYPE kz(As_Array_View(z));
        cg.Solve(
            Make_Krylov_System_Composer< T, KRYLOV_VECTOR_TYPE >(
                Make_Visitor_Sequence(
                    // zero result before GENERIC_SYSTEM_REFERENCE<T>::Apply
                    Make_Compose_Function(FILL1_FUNCTION<T>(0), ARGUMENT_FUNCTION<2>()),
                    PHYSBAM_BOUND_FAST_MEM_FN_TEMPLATE( system, &GENERIC_SYSTEM_REFERENCE<T>::Apply )
                ),
                DOT_PRODUCT_INNER_PRODUCT(),
                MAXIMUM_MAGNITUDE_CONVERGENCE_NORM<T>(),
                NULL_FUNCTION(), // project
                NULL_FUNCTION(), // set_boundary_conditions
                NULL_FUNCTION(), // project_nullspace
                Make_Safe_Jacobi_Preconditioner(
                    PHYSBAM_BOUND_FAST_MEM_FN_TEMPLATE( system, &GENERIC_SYSTEM_REFERENCE<T>::Diag )
                ),
                true, // use_preconditioner
                true  // preconditioner_commutes_with_projection
            ),
            kx, kb, kq, ks, kr, kk, kz,
            tolerance,
            0, // min_iterations
            max_iterations
        );
        lout << "[Executing CONJUGATE_GRADIENT::Solve...] " << timer.Elapsed() << " s" << std::endl;

    }
}

#define EXPLICIT_INSTANTIATION( T ) \
template void Solve_SPD_System_With_ICC_PCG<T>( \
    const unsigned int n_thread, \
    const SOLVER_PARAMS params, \
    const bool has_constant_vectors_in_null_space, \
    const GENERIC_SYSTEM_REFERENCE<T> system, \
    const ARRAY_VIEW<T> rhs, \
    ARRAY_VIEW<T> x, \
    std::ostream& lout);
EXPLICIT_INSTANTIATION( float )
EXPLICIT_INSTANTIATION( double )
#undef EXPLICIT_INSTANTIATION

} // namespace PhysBAM
