//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <cassert>

#include <fstream>
#include <iostream>
#include <string>

#include <boost/preprocessor/seq/enum.hpp>

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_NEGATION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Jeffrey_Utilities/Algorithm/Any_If.h>
#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Eval_Grid_Function.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/CONVERT_INDEX_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/FINE_MULTI_INDEX_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>

#include "Build_Neumann_System.h"
#include "DOMAIN_EMBEDDING_CUBE_SUBSYS.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "DOMAIN_SYSTEM.h"
#include "Evaluate_Error.h"
#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "Params/SOLVER_PARAMS.h"
#include "RAND_MT19937_UNIFORM_REAL.h"

#ifndef PHYSBAM_NO_PETSC
#include <Jeffrey_Utilities/Petsc/CALL_AND_CHKERRQ.h>
#include "Petsc/Solve_SPD_System_With_ICC_PCG.h"
#include "Petsc/Solve_SPD_System_With_ICC_PCG.ipp"
#include <petsc.h>
#endif // #ifndef PHYSBAM_NO_PETSC

#include "Build_And_Solve_Neumann_System.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

namespace
{

template< class T_SYSTEM >
struct INDEX_IS_DIRICHLET;

} // namespace

template< class T, int D >
int Build_And_Solve_Neumann_System(
    const typename EXAMPLE_PARAMS<T,D>::NEUMANN_PARAMS& problem,
    const MAIN_PARAMS<T,D>& main_params,
    typename RAND_MT19937_UNIFORM_REAL<T>::type& rand,
    const ARRAY_VIEW<const T> phi_of_fine_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    typedef DOMAIN_REGULAR_CROSS_SUBSYS<T,D> REGULAR_SUBSYS_TYPE;
    typedef DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D> EMBEDDING_SUBSYS_TYPE;
    typedef DOMAIN_SYSTEM< REGULAR_SUBSYS_TYPE&, EMBEDDING_SUBSYS_TYPE& > SYSTEM_TYPE;

    BASIC_TIMER timer;

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound = As_Vector<int>(main_params.grid.n_cell);
    const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;

    const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
    const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);
    const VECTOR<T,D> dx = (max_x - min_x) / cell_multi_index_bound.max_multi_index;

    std::cout << "Allocating system and rhs...";
    std::cout.flush();
    timer.Restart();
    REGULAR_SUBSYS_TYPE regular_subsys(multi_index_bound, dx);
    EMBEDDING_SUBSYS_TYPE embedding_subsys(multi_index_bound);
    const SYSTEM_TYPE system(regular_subsys, embedding_subsys);
    ARRAY<T> system_rhs(multi_index_bound.Size()); // init'ed to 0
    std::cout << timer.Elapsed() << " s" << std::endl;

    std::cout << "Building Neumann system..." << std::endl;
    timer.Restart();
    Build_Neumann_System(
        problem, main_params,
        phi_of_fine_index,
        regular_subsys, embedding_subsys, ARRAY_VIEW<T>(system_rhs)
    );
    std::cout << "[Building Neumann system...] " << timer.Elapsed() << " s" << std::endl;

    std::cout << "Evaluating residual norm of continuous solution...";
    std::cout.flush();
    timer.Restart();
    ARRAY<T> u_continuous(multi_index_bound.Size()); // init'ed to 0
    Eval_Grid_Function_MT(
        main_params.general.n_thread,
        min_x, max_x, multi_index_bound,
        problem.u,
        ARRAY_VIEW<T>(u_continuous)
    );
    {
        ARRAY<T> residual(-system_rhs);
        system.Apply(u_continuous, residual);
        std::cout << timer.Elapsed() << " s" << std::endl;
        std::cout << "  " << ARRAYS_COMPUTATIONS::Maxabs(residual) << std::endl;
    }

    if(main_params.solver.solver_id == SOLVER_PARAMS::SOLVER_ID_NULL)
        return 0;

    std::cout << "Determining if matrix has a (nontrivial) null space...";
    std::cout.flush();
    timer.Restart();
    const bool has_nontrivial_null_space = !Any_If_MT(
        main_params.general.n_thread,
        1, multi_index_bound.Size(),
        INDEX_IS_DIRICHLET< SYSTEM_TYPE >(system)
    );
    std::cout << timer.Elapsed() << " s" << std::endl;
    std::cout << "  " << (has_nontrivial_null_space ? "yes" : "no") << std::endl;

    std::cout << "Allocating approximate solution...";
    std::cout.flush();
    timer.Restart();
    ARRAY<T> u_approx(multi_index_bound.Size()); // init'ed to 0
    std::cout << timer.Elapsed() << " s" << std::endl;

    switch(main_params.solver.solver_id) {
    case SOLVER_PARAMS::SOLVER_ID_NULL:
        assert(false);
        break;
    case SOLVER_PARAMS::SOLVER_ID_PHYSBAM_CG:
        std::cout << "WARNING: Solver \"physbam-cg\" not yet implemented for Neumann problems." << std::endl;
        break;
    case SOLVER_PARAMS::SOLVER_ID_PHYSBAM_MINRES:
        std::cout << "ERROR: Solver \"physbam-minres\" cannot be used to solve Neumann problems." << std::endl;
        return 1;
#ifdef PHYSBAM_NO_PETSC
    case SOLVER_PARAMS::SOLVER_ID_PETSC_CG:
        std::cout << "WARNING: PETSc not supported on this platform." << std::endl;
        break;
#else // #ifdef PHYSBAM_NO_PETSC
    case SOLVER_PARAMS::SOLVER_ID_PETSC_CG:
        std::cout << "Solving with PETSc CG solver..." << std::endl;
        timer.Restart();
        PHYSBAM_PETSC_CALL_AND_CHKERRQ((
            Petsc::Solve_SPD_System_With_ICC_PCG<T,D>(
                main_params,
                system, system_rhs,
                has_nontrivial_null_space,
                u_approx
            )
        ));
        std::cout << "[Solving with PETSc CG solver...] " << timer.Elapsed() << " s" << std::endl;
        break;
#endif // #ifdef PHYSBAM_NO_PETSC
    case SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES:
        std::cout << "ERROR: Solver \"petsc-minres\" cannot be used to solve Neumann problems." << std::endl;
        return 1;
    case SOLVER_PARAMS::SOLVER_ID_MG:
        std::cout << "WARNING: Solver \"mg\" not yet implemented for Neumann problems." << std::endl;
        break;
    case SOLVER_PARAMS::SOLVER_ID_MGPCG:
        std::cout << "WARNING: Solver \"mgpcg\" not yet implemented for Neumann problems." << std::endl;
        break;
    }

    if(has_nontrivial_null_space) {
        T error_sum = 0;
        int count = 0;
        for(int linear_index = 1; linear_index <= multi_index_bound.Size(); ++linear_index) {
            const MULTI_INDEX_TYPE multi_index = multi_index_bound.Multi_Index(linear_index);
            const MULTI_INDEX_TYPE fine_multi_index = 2 * multi_index - 1;
            const int fine_linear_index = fine_multi_index_bound.Linear_Index(fine_multi_index);
            const T phi = phi_of_fine_index(fine_linear_index);
            const int sign = (0 < phi) - (phi < 0);
            if(sign >= 0)
                continue;
            error_sum += u_approx(linear_index) - u_continuous(linear_index);
            ++count;
        }
        u_approx -= error_sum / count;
    }

    std::cout << "Evaluating error norm in approximate solution...";
    std::cout.flush();
    timer.Restart();
    T max_u_error = 0;
    T max_grad_u_error = 0;
    Evaluate_Error(
        problem, main_params,
        -1,
        Make_Compose_Function(
            SIGN_FUNCTION(),
            phi_of_fine_index,
            fine_multi_index_bound,
            FINE_MULTI_INDEX_FUNCTION<2>(),
            Make_Convert_Index_Function< MULTI_INDEX_TYPE >(multi_index_bound)
        ),
        Make_Compose_Function(
            As_Const_Raw_Array(regular_subsys.sign_of_cell_index),
            cell_multi_index_bound
        ),
        ARRAY_VIEW<const T>(u_approx),
        max_u_error, max_grad_u_error
    );
    std::cout << timer.Elapsed() << " s" << std::endl;
    std::cout << "  |     u_continuous  -      u_approx |_{infty} = " << max_u_error << std::endl;
    std::cout << "  |grad(u_continuous) - grad(u_approx)|_{infty} = " << max_grad_u_error << std::endl;

    {
        const std::string& filename = main_params.output.infty_norm_error_filename;
        if(!filename.empty()) {
            std::ofstream fout(filename.c_str(), std::ios_base::app | std::ios_base::out);
            if(fout.is_open())
                fout << max_u_error << ' ' << max_grad_u_error << std::endl;
            else
                std::cerr << "WARNING: Unable to open file \"" << filename << "\"!" << std::endl;
        }
    }

    return 0;
}

namespace
{

template< class T_SYSTEM >
struct INDEX_IS_DIRICHLET
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        INDEX_IS_DIRICHLET, (( typename T_SYSTEM const &, system ))
    )
public:
    typedef bool result_type;
    bool operator()(const int linear_index) const
    { return system.Stencil_N_Nonzero(linear_index) == 1; }
};

} // namespace

#define EXPLICIT_INSTANTIATION( T, D ) \
template int Build_And_Solve_Neumann_System<T,D>( \
    const EXAMPLE_PARAMS<T,D>::NEUMANN_PARAMS& problem, \
    const MAIN_PARAMS<T,D>& main_params, \
    RAND_MT19937_UNIFORM_REAL<T>::type& rand, \
    const ARRAY_VIEW<const T> phi_of_fine_index);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM
