//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <cassert>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include <boost/function.hpp>
#include <boost/preprocessor/seq/enum.hpp>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/Eval_Grid_Function.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>
#include <Jeffrey_Utilities/GENERIC_SYSTEM_REFERENCE.h>
#include <Jeffrey_Utilities/Grid/ASSIGN_SIGN_TO_INDEX_GRID_VISITOR.h>
#include <Jeffrey_Utilities/Grid/VERTEX_IS_COARSE_CELL_CENTER.h>
#include <Jeffrey_Utilities/Grid/Visit_Cells_With_Sign_Via_Fine_Vertex_Sign.h>
#include <Jeffrey_Utilities/Has_Constant_Vectors_In_Null_Space.h>
#include <Jeffrey_Utilities/Krylov/Solve_SPD_System_With_ICC_PCG.h>
#include <Jeffrey_Utilities/Multi_Index/FINE_MULTI_INDEX_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Stencils/INDEX_TRANSFORM_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_DIFFERENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_NEGATION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "AGGREGATE_CONSTRAINT_SYSTEM.h"
#include "Aggregate_Constraints.h"
#include "Build_Dirichlet_System.h"
#include "Build_ZTAZ_Embedding_Subsys.h"
#include "DIRICHLET_CONSTRAINT_SYSTEM.h"
#include "DOMAIN_EMBEDDING_CUBE_SUBSYS.h"
#include "DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "DOMAIN_SYSTEM.h"
#include "Evaluate_Error.h"
#include "Init_ZTAZ_Embedding.h"
#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "RAND_MT19937_UNIFORM_REAL.h"
#include "Select_Indys.h"

#ifdef PHYSBAM_USE_PETSC
#include <petsc.h>
#include <Jeffrey_Utilities/Petsc/CALL_AND_CHKERRQ.h>
#include <Jeffrey_Utilities/Petsc/Solve_SPD_System_With_ICC_PCG.h>
#endif // #ifdef PHYSBAM_USE_PETSC

#include "Build_And_Solve_Dirichlet_System.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
int Build_And_Solve_Dirichlet_System(
    typename EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    typename RAND_MT19937_UNIFORM_REAL<T>::type& rand,
    const ARRAY_VIEW<const T> phi_of_fine_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    typedef DOMAIN_REGULAR_CROSS_SUBSYS<T,D> REGULAR_SUBSYS_TYPE;
    typedef DOMAIN_EMBEDDING_CUBE_SUBSYS<T,D> EMBEDDING_SUBSYS_TYPE;
    typedef DOMAIN_SYSTEM< REGULAR_SUBSYS_TYPE&, EMBEDDING_SUBSYS_TYPE& > SYSTEM_TYPE;
    typedef DIRICHLET_CONSTRAINT_SYSTEM<T,D> CONSTRAINT_SYSTEM_TYPE;

    typedef DOMAIN_EMBEDDING_UNSTRUCTURED_SUBSYS<T,D> ZTAZ_EMBEDDING_SUBSYS_TYPE;
    typedef DOMAIN_SYSTEM< REGULAR_SUBSYS_TYPE&, ZTAZ_EMBEDDING_SUBSYS_TYPE& > ZTAZ_SYSTEM_TYPE;
    typedef AGGREGATE_CONSTRAINT_SYSTEM<T,D> AGGREGATE_CONSTRAINT_SYSTEM_TYPE;

    BASIC_TIMER timer;

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound = As_Vector<int>(main_params.grid.n_cell);
    const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;

    const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
    const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);
    const VECTOR<T,D> dx = (max_x - min_x) / cell_multi_index_bound.max_multi_index;

    typename Result_Of::MAKE_COMPOSE_FUNCTION<
        SIGN_FUNCTION, ARRAY_VIEW<const T>, MULTI_INDEX_BOUND<D>
    >::type const sign_of_fine_index = Make_Compose_Function(
        SIGN_FUNCTION(), phi_of_fine_index, fine_multi_index_bound
    );

    std::cout << "Allocating system and rhs...";
    std::cout.flush();
    timer.Restart();
    REGULAR_SUBSYS_TYPE regular_subsys(multi_index_bound, dx);
    EMBEDDING_SUBSYS_TYPE embedding_subsys(multi_index_bound);
    const SYSTEM_TYPE system(regular_subsys, embedding_subsys);
    CONSTRAINT_SYSTEM_TYPE constraint_system(multi_index_bound);
    ARRAY<T> system_rhs(multi_index_bound.Size()); // init'ed to 0
    ARRAY<T> constraint_rhs;
    std::cout << timer.Elapsed() << " s" << std::endl;

    std::cout << "Initializing sign_of_cell_index...";
    std::cout.flush();
    ARRAY<signed char> sign_of_cell_index(cell_multi_index_bound.Size()); // init'ed to 0
    Visit_Cells_With_Sign_Via_Fine_Vertex_Sign_MT<2>(
        main_params.general.n_thread,
        cell_multi_index_bound,
        sign_of_fine_index,
        Make_Assign_Sign_To_Index_Grid_Visitor(
            Make_Array_Wrapper_Function(sign_of_cell_index)
        ),
        -1 // sign_of_zero
    );
    std::cout << timer.Elapsed() << " s" << std::endl;

    std::cout << "Building Dirichlet system..." << std::endl;
    timer.Restart();
    Build_Dirichlet_System(
        problem, main_params,
        phi_of_fine_index,
        As_Const_Array_View(sign_of_cell_index),
        regular_subsys, embedding_subsys, As_Array_View(system_rhs),
        constraint_system, constraint_rhs,
        std::cout
    );
    std::cout << "[Building Dirichlet system...] " << timer.Elapsed() << " s" << std::endl;

    const int n_embedding = embedding_subsys.stencils.Size();
    const int n_constraint = constraint_system.stencils.Size();

    std::cout << "Evaluating constraint residual norm of continuous solution...";
    std::cout.flush();
    timer.Restart();
    ARRAY<T> u_continuous(multi_index_bound.Size()); // init'ed to 0
    Eval_Grid_Function_MT(
        main_params.general.n_thread,
        min_x, max_x, multi_index_bound,
        problem.u,
        As_Array_View(u_continuous)
    );
    {
        ARRAY<T> constraint_residual(-constraint_rhs);
        constraint_system.Apply(u_continuous, constraint_residual);
        std::cout << timer.Elapsed() << " s" << std::endl;
        std::cout << "  " << ARRAYS_COMPUTATIONS::Maxabs(constraint_residual) << std::endl;
    }

    if(
        main_params.solver.solver_id == SOLVER_PARAMS::SOLVER_ID_NULL
     && problem.constraint_id == EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_NULL
    )
        return 0;

    std::cout << "Allocating approximate solution...";
    std::cout.flush();
    timer.Restart();
    ARRAY<T> u_approx(multi_index_bound.Size()); // init'ed to 0
    std::cout << timer.Elapsed() << " s" << std::endl;

    if(problem.constraint_id == EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_SINGLE_CELL) {
        switch(main_params.solver.solver_id) {
        case SOLVER_PARAMS::SOLVER_ID_PHYSBAM_MINRES:
            std::cout << "WARNING: Solver \"physbam-minres\" not yet implemented for Dirichlet problems." << std::endl;
            break;
        case SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES:
#ifdef PHYSBAM_USE_PETSC
            std::cout << "WARNING: Solver \"petsc-minres\" not yet implemented for Dirichlet problems." << std::endl;
#else // #ifdef PHYSBAM_USE_PETSC
            std::cout << "WARNING: PETSc not supported on this platform." << std::endl;
#endif // #ifdef PHYSBAM_USE_PETSC
            break;
        default:
            std::cout << "ERROR: Must use either \"physbam-minres\" or \"petsc-minres\" "
                         "solvers for Dirichlet problems with single-cell constraints."
                      << std::endl;
            return 1;
        }
    }
    else {

        std::cout << "Selecting indys...";
        std::cout.flush();
        timer.Restart();
        ARRAY<int> indy_index_of_constraint_index(n_constraint, false); // uninit'ed
        AGGREGATE_CONSTRAINT_SYSTEM_TYPE aggregate_constraint_system;
        aggregate_constraint_system.index_of_indy_index.Preallocate(n_embedding / (1 << D));
        {
            boost::function< bool ( int ) > indyable;
            switch(problem.constraint_id) {
            case EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_DOUBLE_CELL:
                indyable = Make_Vertex_Is_Coarse_Cell_Center<2>(multi_index_bound);
                break;
            case EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_AGGREGATE:
                indyable = Make_Compose_Function(
                    Make_Equal_Function(+1),
                    sign_of_fine_index,
                    FINE_MULTI_INDEX_FUNCTION<2>(),
                    multi_index_bound
                );
                break;
            default:
                assert(false);
            }
            Select_Indys(
                cell_multi_index_bound,
                multi_index_bound,
                constraint_system.stencil_index_of_cell_linear_index,
                Make_Compose_Function(
                    Make_Index_Transform_Stencil_Proxy_Function(multi_index_bound),
                    BOUND_FAST_MEM_FN<
                        typename CONSTRAINT_SYSTEM_TYPE::CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
                            (CONSTRAINT_SYSTEM_TYPE::*)( int ) const,
                        &CONSTRAINT_SYSTEM_TYPE::Multi_Index_Stencil_Proxy
                    >(constraint_system)
                ),
                indyable,
                problem.min_relative_indy_weight,
                As_Array_View(indy_index_of_constraint_index),
                aggregate_constraint_system.index_of_indy_index
            );
        }
        std::cout << timer.Elapsed() << " s" << std::endl;
        const int n_indy = aggregate_constraint_system.index_of_indy_index.Size();
        std::cout << "  # of indys = " << n_indy << std::endl;

        std::cout << "Aggregating constraints...";
        std::cout.flush();
        timer.Restart();
        aggregate_constraint_system.Init_Indy_Index_Of_Index();
        aggregate_constraint_system.value_of_indy_index.Exact_Resize(n_indy); // init'ed to 0
        aggregate_constraint_system.indyless_stencils.Exact_Resize(n_indy);
        Aggregate_Constraints(
            Make_Compose_Function(
                Make_Index_Transform_Stencil_Proxy_Function(multi_index_bound),
                BOUND_FAST_MEM_FN<
                    typename CONSTRAINT_SYSTEM_TYPE::CONST_MULTI_INDEX_STENCIL_PROXY_TYPE
                        (CONSTRAINT_SYSTEM_TYPE::*)( int ) const,
                    &CONSTRAINT_SYSTEM_TYPE::Multi_Index_Stencil_Proxy
                >(constraint_system)
            ),
            As_Const_Array_View(indy_index_of_constraint_index),
            As_Const_Array_View(aggregate_constraint_system.index_of_indy_index),
            As_Array_View(aggregate_constraint_system.value_of_indy_index),
            BOUND_FAST_MEM_FN<
                typename AGGREGATE_CONSTRAINT_SYSTEM_TYPE::INDYLESS_STENCIL_PROXY_TYPE
                    (AGGREGATE_CONSTRAINT_SYSTEM_TYPE::*)( int ),
                &AGGREGATE_CONSTRAINT_SYSTEM_TYPE::Indyless_Stencil_Proxy
            >(aggregate_constraint_system)
        );
        aggregate_constraint_system.Init_Stencils_Containing_Index();
        ARRAY<T> aggregate_constraint_rhs(n_indy); // init'ed to 0
        for(int constraint_index = 1; constraint_index <= n_constraint; ++constraint_index) {
            const int indy_index = indy_index_of_constraint_index(constraint_index);
            aggregate_constraint_rhs(indy_index) += constraint_rhs(constraint_index);
        }
        std::cout << timer.Elapsed() << " s" << std::endl;
        if(n_indy != 0) {
            const T max_indy_weight = aggregate_constraint_system.value_of_indy_index(1);
            const T min_indy_weight = aggregate_constraint_system.value_of_indy_index(n_indy);
            std::cout << "  max indy weight = " << max_indy_weight << '\n'
                      << "  min indy weight = " << min_indy_weight << '\n'
                      << "    ratio = " << max_indy_weight / min_indy_weight
                      << std::endl;
        }

        std::cout << "Evaluating aggregate constraint residual norm of continuous solution...";
        std::cout.flush();
        timer.Restart();
        {
            ARRAY<T> aggregate_constraint_residual(-aggregate_constraint_rhs);
            aggregate_constraint_system.Apply(u_continuous, aggregate_constraint_residual);
            std::cout << timer.Elapsed() << " s" << std::endl;
            std::cout << "  " << ARRAYS_COMPUTATIONS::Maxabs(aggregate_constraint_residual) << std::endl;
        }

        if(
            main_params.solver.solver_id == SOLVER_PARAMS::SOLVER_ID_PHYSBAM_MINRES
         || main_params.solver.solver_id == SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES
        ) {
            switch(main_params.solver.solver_id) {
            case SOLVER_PARAMS::SOLVER_ID_PHYSBAM_MINRES:
                std::cout << "WARNING: Solver \"physbam-minres\" not yet implemented for Dirichlet problems." << std::endl;
                break;
            case SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES:
#ifdef  PHYSBAM_USE_PETSC
                std::cout << "WARNING: Solver \"petsc-minres\" not yet implemented for Dirichlet problems." << std::endl;
#else // #ifdef  PHYSBAM_USE_PETSC
                std::cout << "WARNING: PETSc not supported on this platform." << std::endl;
                break;
#endif // #ifdef  PHYSBAM_USE_PETSC
            default:
                assert(false);
            }
        }
        else {

            std::cout << "Building Z^T*A*Z embedding subsystem...";
            std::cout.flush();
            timer.Restart();
            ZTAZ_EMBEDDING_SUBSYS_TYPE ztaz_embedding_subsys(multi_index_bound);
            Init_ZTAZ_Embedding(
                embedding_subsys.stencil_index_of_linear_index,
                As_Const_Array_View(embedding_subsys.linear_index_of_stencil_index),
                aggregate_constraint_system.indy_index_of_index,
                As_Const_Array_View(aggregate_constraint_system.index_of_indy_index),
                Make_Compose_Function(
                    BOUND_FAST_MEM_FN<
                        typename REGULAR_SUBSYS_TYPE::CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
                            (REGULAR_SUBSYS_TYPE::*)( int ) const,
                        &REGULAR_SUBSYS_TYPE::Linear_Index_Stencil_Proxy
                    >(regular_subsys),
                    As_Const_Array_View(aggregate_constraint_system.index_of_indy_index)
                ),
                ztaz_embedding_subsys.stencil_index_of_linear_index,
                ztaz_embedding_subsys.linear_index_of_stencil_index
            );
            const int n_ztaz_embedding = ztaz_embedding_subsys.stencil_index_of_linear_index.Size();
            ztaz_embedding_subsys.stencils.Exact_Resize(n_ztaz_embedding);
            Build_ZTAZ_Embedding_Subsys(
                BOUND_FAST_MEM_FN<
                    typename REGULAR_SUBSYS_TYPE::CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
                        (REGULAR_SUBSYS_TYPE::*)( int ) const,
                    &REGULAR_SUBSYS_TYPE::Linear_Index_Stencil_Proxy
                >(regular_subsys),
                BOUND_FAST_MEM_FN<
                    typename EMBEDDING_SUBSYS_TYPE::CONST_LINEAR_INDEX_STENCIL_PROXY_TYPE
                        (EMBEDDING_SUBSYS_TYPE::*)( int ) const,
                    &EMBEDDING_SUBSYS_TYPE::Linear_Index_Stencil_Proxy
                >(embedding_subsys),
                aggregate_constraint_system.indy_index_of_index,
                As_Const_Array_View(aggregate_constraint_system.index_of_indy_index),
                aggregate_constraint_system.stencils_containing_index,
                As_Const_Array_View(aggregate_constraint_system.value_of_indy_index),
                BOUND_FAST_MEM_FN<
                    typename AGGREGATE_CONSTRAINT_SYSTEM_TYPE::CONST_INDYLESS_STENCIL_PROXY_TYPE
                        (AGGREGATE_CONSTRAINT_SYSTEM_TYPE::*)( int ) const,
                    &AGGREGATE_CONSTRAINT_SYSTEM_TYPE::Indyless_Stencil_Proxy
                >(aggregate_constraint_system),
                As_Const_Array_View(ztaz_embedding_subsys.linear_index_of_stencil_index),
                BOUND_FAST_MEM_FN<
                    typename ZTAZ_EMBEDDING_SUBSYS_TYPE::LINEAR_INDEX_STENCIL_PROXY_TYPE
                        (ZTAZ_EMBEDDING_SUBSYS_TYPE::*)( int ),
                    &ZTAZ_EMBEDDING_SUBSYS_TYPE::Linear_Index_Stencil_Proxy_Of_Stencil_Index
                >(ztaz_embedding_subsys)
            );
            std::cout << timer.Elapsed() << " s" << std::endl;

#if 0
            {
                std::cout << "Verifying correctness of Z^T*A*Z system via multiplication by standard basis vectors...";
                std::cout.flush();
                timer.Restart();
                ARRAY<T> e(multi_index_bound.Size()); // init'ed to 0
                ARRAY<T> ztaz_e_a(multi_index_bound.Size(), false); // uninit'ed
                ARRAY<T> ztaz_e_b(multi_index_bound.Size(), false); // uninit'ed
                T max_abs_diff = 0;
                for(int linear_index = 1; linear_index <= multi_index_bound.Size(); ++linear_index) {
                    if(aggregate_constraint_system.indy_index_of_index.Contains(linear_index))
                        continue;
                    e(linear_index) = static_cast<T>(1);
                    {
                        ARRAY<T>& z_e = ztaz_e_b;
                        z_e = e;
                        aggregate_constraint_system.Apply_Z(z_e);
                        ztaz_e_a.Fill(static_cast<T>(0));
                        system.Apply(z_e, ztaz_e_a);
                        aggregate_constraint_system.Apply_Z_Transpose(ztaz_e_a);
                    }
                    ztaz_e_b.Fill(static_cast<T>(0));
                    regular_subsys.Apply(e, ztaz_e_b);
                    ztaz_embedding_subsys.Apply(e, ztaz_e_b);
                    for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
                        const int indy_linear_index = aggregate_constraint_system.index_of_indy_index(indy_index);
                        ztaz_e_b(indy_linear_index) = static_cast<T>(0);
                    }
                    const T local_max_abs_diff = ARRAYS_COMPUTATIONS::Maxabs(ztaz_e_a - ztaz_e_b);
                    max_abs_diff = std::max(max_abs_diff, local_max_abs_diff);
                    e(linear_index) = static_cast<T>(0);
                }
                std::cout << timer.Elapsed() << " s" << std::endl;
                std::cout << "  max abs diff = " << max_abs_diff << std::endl;
            }
#endif // #if 0|1

            if(main_params.general.randomized_check) {
                std::cout << "Verifying correctness of Z^T*A*Z system via multiplication by a random vector...";
                std::cout.flush();
                timer.Restart();
                ARRAY<T> u_random(multi_index_bound.Size()); // init'ed to 0
                for(int linear_index = 1; linear_index <= multi_index_bound.Size(); ++linear_index)
                    if(
                        system.Stencil_N_Nonzero(linear_index) != 0
                     && !aggregate_constraint_system.indy_index_of_index.Contains(linear_index)
                    )
                        u_random(linear_index) = rand();
                ARRAY<T> ztaz_u_random_a(multi_index_bound.Size()); // init'ed to 0
                {
                    ARRAY<T> z_u_random(u_random);
                    aggregate_constraint_system.Apply_Z(z_u_random);
                    system.Apply(z_u_random, ztaz_u_random_a);
                    aggregate_constraint_system.Apply_Z_Transpose(ztaz_u_random_a);
                }
                ARRAY<T> ztaz_u_random_b(multi_index_bound.Size()); // init'ed to 0
                regular_subsys.Apply(u_random, ztaz_u_random_b);
                ztaz_embedding_subsys.Apply(u_random, ztaz_u_random_b);
                for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
                    const int indy_linear_index = aggregate_constraint_system.index_of_indy_index(indy_index);
                    ztaz_u_random_b(indy_linear_index) = static_cast<T>(0);
                }
                const T max_abs_diff = ARRAYS_COMPUTATIONS::Maxabs(ztaz_u_random_a - ztaz_u_random_b);
                std::cout << timer.Elapsed() << " s" << std::endl;
                std::cout << "  max abs diff = " << max_abs_diff << std::endl;
            }

            std::cout << "Evaluating Z^T*(f - A*c)...";
            std::cout.flush();
            timer.Restart();
            ARRAY<T>& ztaz_system_rhs = system_rhs;
            for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
                const int linear_index = aggregate_constraint_system.index_of_indy_index(indy_index);
                const T c_value = aggregate_constraint_rhs(indy_index)
                                / aggregate_constraint_system.value_of_indy_index(indy_index);
                system.Apply_Transpose(linear_index, ztaz_system_rhs, -c_value);
            }
            aggregate_constraint_system.Apply_Z_Transpose(ztaz_system_rhs);
            std::cout << timer.Elapsed() << " s" << std::endl;

            std::cout << "Zero'ing indy stencils in Z^T*A*Z regular subsystem...";
            std::cout.flush();
            timer.Restart();
            for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
                const int linear_index = aggregate_constraint_system.index_of_indy_index(indy_index);
                assert(!ztaz_embedding_subsys.stencil_index_of_linear_index.Contains(linear_index));
                regular_subsys.Zero_Stencil(linear_index);
            }
            std::cout << timer.Elapsed() << " s" << std::endl;

            const ZTAZ_SYSTEM_TYPE ztaz_system(regular_subsys, ztaz_embedding_subsys);

            std::cout << "Evaluating Z^T*A*Z residual norm of continuous solution...";
            std::cout.flush();
            timer.Restart();
            {
                ARRAY<T> ztaz_residual(-ztaz_system_rhs);
                ztaz_system.Apply(u_continuous, ztaz_residual);
                std::cout << timer.Elapsed() << " s" << std::endl;
                std::cout << "  " << ARRAYS_COMPUTATIONS::Maxabs(ztaz_residual) << std::endl;
            }

            std::cout << "Determining if constant vectors in null space...";
            std::cout.flush();
            timer.Restart();
            const bool has_constant_vectors_in_null_space =
                Has_Constant_Vectors_In_Null_Space<T>(
                    main_params.general.n_thread,
                    multi_index_bound.Size(),
                    ztaz_system
                );
            std::cout << timer.Elapsed() << " s" << std::endl;
            std::cout << "  " << (has_constant_vectors_in_null_space ? "yes" : "no") << std::endl;

            switch(main_params.solver.solver_id) {
            case SOLVER_PARAMS::SOLVER_ID_NULL:
                return 0;
            case SOLVER_PARAMS::SOLVER_ID_PHYSBAM_CG:
                std::cout << "Solving with PhysBAM CG solver..." << std::endl;
                timer.Restart();
                PhysBAM::Solve_SPD_System_With_ICC_PCG(
                    main_params.general.n_thread,
                    main_params.solver,
                    has_constant_vectors_in_null_space,
                    GENERIC_SYSTEM_REFERENCE<T>(ztaz_system),
                    As_Array_View(ztaz_system_rhs),
                    As_Array_View(u_approx),
                    std::cout
                );
                std::cout << "[Solving with PhysBAM CG solver...] " << timer.Elapsed() << " s" << std::endl;
                break;
            case SOLVER_PARAMS::SOLVER_ID_PETSC_CG:
#ifdef PHYSBAM_USE_PETSC
                std::cout << "Solving with PETSc CG solver..." << std::endl;
                timer.Restart();
                PHYSBAM_PETSC_CALL_AND_CHKERRQ((
                    Petsc::Solve_SPD_System_With_ICC_PCG(
                        main_params.general.n_thread,
                        main_params.solver,
                        has_constant_vectors_in_null_space,
                        GENERIC_SYSTEM_REFERENCE<T>(ztaz_system),
                        As_Const_Array_View(ztaz_system_rhs),
                        As_Array_View(u_approx),
                        std::cout
                    )
                ));
                std::cout << "[Solving with PETSc CG solver...] " << timer.Elapsed() << " s" << std::endl;
#else // #ifdef PHYSBAM_USE_PETSC
                std::cout << "WARNING: PETSc not supported on this platform." << std::endl;
#endif // #ifdef PHYSBAM_USE_PETSC
                break;
            case SOLVER_PARAMS::SOLVER_ID_MG:
                std::cout << "WARNING: Solver \"mg\" not yet implemented for Dirichlet problems." << std::endl;
                break;
            case SOLVER_PARAMS::SOLVER_ID_MGPCG:
                std::cout << "WARNING: Solver \"mgpcg\" not yet implemented for Dirichlet problems." << std::endl;
                break;
            default:
                assert(false);
            }

            std::cout << "Evaluating c + Z*v...";
            std::cout.flush();
            timer.Restart();
            aggregate_constraint_system.Apply_Z(u_approx);
            for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
                const int linear_index = aggregate_constraint_system.index_of_indy_index(indy_index);
                const T c_value = aggregate_constraint_rhs(indy_index)
                                / aggregate_constraint_system.value_of_indy_index(indy_index);
                u_approx(linear_index) += c_value;
            }
            std::cout << timer.Elapsed() << " s" << std::endl;

        }
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
            sign_of_fine_index,
            FINE_MULTI_INDEX_FUNCTION<2>(),
            multi_index_bound
        ),
        Make_Compose_Function(
            As_Const_Array_View(sign_of_cell_index),
            cell_multi_index_bound
        ),
        As_Const_Array_View(u_approx),
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

#define EXPLICIT_INSTANTIATION( T, D ) \
template int Build_And_Solve_Dirichlet_System<T,D>( \
    EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS const & problem, \
    const MAIN_PARAMS<T,D>& main_params, \
    RAND_MT19937_UNIFORM_REAL<T>::type& rand, \
    const ARRAY_VIEW<const T> phi_of_fine_index);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
