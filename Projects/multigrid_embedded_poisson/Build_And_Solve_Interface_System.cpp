//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <cmath>

#include <fstream>
#include <iostream>

#include <boost/function.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/mpl/vector/vector10.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/ref.hpp>

#include <Jeffrey_Utilities/Algorithm/For_Each.h>
#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/Functional/APPLY_AND_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/APPLY_ASSIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/CONSTANT_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/IF_ELSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Grid/VERTEX_IS_COARSE_CELL_CENTER.h>
#include <Jeffrey_Utilities/Multi_Index/FINE_MULTI_INDEX_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_X_FUNCTION.h>
#include <Jeffrey_Utilities/Stencils/INDEX_TRANSFORM_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/Stencils/ZERO_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_DIFFERENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_NEGATION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "AGGREGATE_CONSTRAINT_SYSTEM.h"
#include "Aggregate_Constraints.h"
#include "Build_Interface_System.h"
#include "Build_ZTAZ_Embedding_Subsys.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "EMBEDDING_UNSTRUCTURED_SUBSYS.h"
#include "Evaluate_Error.h"
#include "Init_ZTAZ_Embedding.h"
#include "INTERFACE_CONSTRAINT_SYSTEM.h"
#include "INTERFACE_INDEX_TRANSFORM.h"
#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "RAND_MT19937_UNIFORM_REAL.h"
#include "Select_Indys.h"
#include "SYSTEM_SUM.h"

#ifndef PHYSBAM_NO_PETSC
#include <Jeffrey_Utilities/Petsc/CALL_AND_CHKERRQ.h>
#include "Petsc/Solve_SPD_System_With_ICC_PCG.h"
#include "Petsc/Solve_SPD_System_With_ICC_PCG.ipp"
#include <petsc.h>
#endif // #ifndef PHYSBAM_NO_PETSC

#include "Build_And_Solve_Interface_System.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
int Build_And_Solve_Interface_System(
    const typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS& problem,
    const MAIN_PARAMS<T,D>& main_params,
    typename RAND_MT19937_UNIFORM_REAL<T>::type& rand,
    const ARRAY_VIEW<const T> phi_of_fine_index)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    typedef DOMAIN_REGULAR_CROSS_SUBSYS<T,D> REGULAR_SUBSYS_TYPE;
    typedef EMBEDDING_UNSTRUCTURED_SUBSYS<T> EMBEDDING_SUBSYS_TYPE;
    typedef SYSTEM_SUM< boost::mpl::vector2< REGULAR_SUBSYS_TYPE&, EMBEDDING_SUBSYS_TYPE& > > SYSTEM_TYPE;
    typedef INTERFACE_CONSTRAINT_SYSTEM<T,D> CONSTRAINT_SYSTEM_TYPE;

    typedef EMBEDDING_UNSTRUCTURED_SUBSYS<T> ZTAZ_EMBEDDING_SUBSYS_TYPE;
    typedef SYSTEM_SUM< boost::mpl::vector2< REGULAR_SUBSYS_TYPE&, ZTAZ_EMBEDDING_SUBSYS_TYPE& > > ZTAZ_SYSTEM_TYPE;
    typedef AGGREGATE_CONSTRAINT_SYSTEM<T,D> AGGREGATE_CONSTRAINT_SYSTEM_TYPE;

    BASIC_TIMER timer;

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound = As_Vector<int>(main_params.grid.n_cell);
    const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;

    const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
    const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);
    const VECTOR<T,D> dx = (max_x - min_x) / cell_multi_index_bound.max_multi_index;

    const MULTI_INDEX_X_FUNCTION< T, D, MULTI_INDEX_BOUND<D> >
        x_of_index(min_x, max_x, multi_index_bound);

    std::cout << "Allocating system...";
    std::cout.flush();
    timer.Restart();
    REGULAR_SUBSYS_TYPE regular_subsys(multi_index_bound, dx);
    EMBEDDING_SUBSYS_TYPE embedding_subsys;
    SYSTEM_TYPE system(boost::fusion::make_vector(
        boost::ref(regular_subsys),
        boost::ref(embedding_subsys))
    );
    CONSTRAINT_SYSTEM_TYPE constraint_system;
    ARRAY<T> system_rhs;
    ARRAY<T> constraint_rhs;
    std::cout << timer.Elapsed() << " s" << std::endl;

    std::cout << "Building interface system..." << std::endl;
    timer.Restart();
    Build_Interface_System(
        problem, main_params,
        phi_of_fine_index,
        regular_subsys, embedding_subsys, system_rhs,
        constraint_system, constraint_rhs,
        std::cout
    );
    std::cout << "[Building interface system...] " << timer.Elapsed() << " s" << std::endl;

    const int n_embedding  = embedding_subsys.stencils.Size();
    const int n_virtual    = n_embedding / 2;
    const int n_index      = multi_index_bound.Size() + n_virtual;
    const int n_constraint = constraint_system.stencils.Size();

    // Construct index_transform.
    typedef BOUND_FAST_MEM_FN<
        const int& (HASHTABLE<int,int>::*)( const int& ) const,
        &HASHTABLE<int,int>::Get
    > VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX_TYPE;
    typedef ARRAY_VIEW<const int> GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET_TYPE;
    typedef typename Result_Of::MAKE_COMPOSE_FUNCTION<
        SIGN_FUNCTION,
        ARRAY_VIEW<const T>,
        MULTI_INDEX_BOUND<D>,
        FINE_MULTI_INDEX_FUNCTION<2>,
        MULTI_INDEX_BOUND<D>
    >::type SIGN_OF_GRID_INDEX_TYPE;
    typedef INTERFACE_INDEX_TRANSFORM<
        VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX_TYPE,
        GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET_TYPE,
        SIGN_OF_GRID_INDEX_TYPE
    > INDEX_TRANSFORM_TYPE;
    const INDEX_TRANSFORM_TYPE index_transform(
        multi_index_bound.Size(),
        VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX_TYPE(embedding_subsys.stencil_index_of_index),
        GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET_TYPE(embedding_subsys.index_of_stencil_index),
        Make_Compose_Function(
            SIGN_FUNCTION(),
            phi_of_fine_index,
            fine_multi_index_bound,
            FINE_MULTI_INDEX_FUNCTION<2>(),
            multi_index_bound
        )
    );

    BOUND_FAST_MEM_FN<
        int (INDEX_TRANSFORM_TYPE::*)( int ) const,
        &INDEX_TRANSFORM_TYPE::Grid_Index_Of_Index
    > grid_index_of_index(index_transform);

    std::cout << "Evaluating constraint residual norm of continuous solution...";
    std::cout.flush();
    timer.Restart();
    ARRAY<T> u_continuous(multi_index_bound.Size() + n_embedding/2); // init'ed to 0
    For_Each_MT(
        main_params.general.n_thread,
        1, n_index,
        Make_Apply_Assign_Function(
            Make_Array_Wrapper_Function(u_continuous),
            Make_If_Else_Function(
                Make_Compose_Function(
                    Make_Equal_Function(-1),
                    PHYSBAM_BOUND_FAST_MEM_FN_TEMPLATE(
                        index_transform,
                        &INDEX_TRANSFORM_TYPE::Domain_Sign_Of_Index
                    )
                ),
                Make_Compose_Function(problem.negative.u, x_of_index, grid_index_of_index),
                Make_Compose_Function(problem.positive.u, x_of_index, grid_index_of_index)
            )
        )
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
    ARRAY<T> u_approx(n_index); // init'ed to 0
    std::cout << timer.Elapsed() << " s" << std::endl;

    if(problem.constraint_id == EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_SINGLE_CELL) {
        switch(main_params.solver.solver_id) {
        case SOLVER_PARAMS::SOLVER_ID_PHYSBAM_MINRES:
            std::cout << "WARNING: Solver \"physbam-minres\" not yet implemented for interface problems." << std::endl;
            break;
#ifdef PHYSBAM_NO_PETSC
        case SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES:
            std::cout << "WARNING: PETSc not supported on this platform." << std::endl;
            break;
#else // #ifdef PHYSBAM_NO_PETSC
        case SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES:
            std::cout << "WARNING: Solver \"petsc-minres\" not yet implemented for interface problems." << std::endl;
            break;
#endif // #ifdef PHYSBAM_NO_PETSC
        default:
            std::cout << "ERROR: Must use either \"physbam-minres\" or \"petsc-minres\" "
                         "solvers for interface problems with single-cell constraints."
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
        aggregate_constraint_system.index_of_indy_index.Preallocate(n_embedding / (2 * (1 << D)));
        {
            boost::function< bool ( int ) > indyable;
            const BOUND_FAST_MEM_FN<
                bool (INDEX_TRANSFORM_TYPE::*)( int ) const,
                &INDEX_TRANSFORM_TYPE::Index_Is_Virtual
            > index_is_virtual(index_transform);
            switch(problem.constraint_id) {
            case EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_DOUBLE_CELL:
                indyable = Make_Apply_And_Function(
                    index_is_virtual,
                    Make_Compose_Function(
                        Make_Vertex_Is_Coarse_Cell_Center<2>(multi_index_bound),
                        grid_index_of_index
                    )
                );
                break;
            case EXAMPLE_PARAMS_BASE::CONSTRAINT_ID_AGGREGATE:
                indyable = index_is_virtual;
                break;
            default:
                assert(false);
            }
            Select_Indys(
                cell_multi_index_bound,
                Make_Compose_Function(multi_index_bound, grid_index_of_index),
                constraint_system.stencil_index_of_cell_index,
                BOUND_FAST_MEM_FN<
                    typename CONSTRAINT_SYSTEM_TYPE::CONST_STENCIL_PROXY_TYPE
                        (CONSTRAINT_SYSTEM_TYPE::*)( int ) const,
                    &CONSTRAINT_SYSTEM_TYPE::Stencil_Proxy
                >(constraint_system),
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
            BOUND_FAST_MEM_FN<
                typename CONSTRAINT_SYSTEM_TYPE::CONST_STENCIL_PROXY_TYPE
                    (CONSTRAINT_SYSTEM_TYPE::*)( int ) const,
                &CONSTRAINT_SYSTEM_TYPE::Stencil_Proxy
            >(constraint_system),
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
            const T max_indy_weight = std::abs(aggregate_constraint_system.value_of_indy_index(1));
            const T min_indy_weight = std::abs(aggregate_constraint_system.value_of_indy_index(n_indy));
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
#ifdef PHYSBAM_NO_PETSC
            case SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES:
                std::cout << "WARNING: PETSc not supported on this platform." << std::endl;
                break;
#else // #ifdef PHYSBAM_NO_PETSC
            case SOLVER_PARAMS::SOLVER_ID_PETSC_MINRES:
                std::cout << "WARNING: Solver \"petsc-minres\" not yet implemented for Dirichlet problems." << std::endl;
                break;
#endif // #ifdef PHYSBAM_NO_PETSC
            default:
                assert(false);
            }
        }
        else {

            std::cout << "Building Z^T*A*Z embedding subsystem...";
            std::cout.flush();
            timer.Restart();
            ZTAZ_EMBEDDING_SUBSYS_TYPE ztaz_embedding_subsys;
            Init_ZTAZ_Embedding(
                embedding_subsys.stencil_index_of_index,
                As_Const_Array_View(embedding_subsys.index_of_stencil_index),
                aggregate_constraint_system.indy_index_of_index,
                As_Const_Array_View(aggregate_constraint_system.index_of_indy_index),
                Make_Constant_Function(ZERO_STENCIL_PROXY<int,T>()), // okay, since all indys are virtual
                ztaz_embedding_subsys.stencil_index_of_index,
                ztaz_embedding_subsys.index_of_stencil_index
            );
            const int n_ztaz_embedding = ztaz_embedding_subsys.stencil_index_of_index.Size();
            ztaz_embedding_subsys.stencils.Exact_Resize(n_ztaz_embedding);
            Build_ZTAZ_Embedding_Subsys(
                Make_Constant_Function(ZERO_STENCIL_PROXY<int,T>()), // okay, since all indys are virtual
                BOUND_FAST_MEM_FN<
                    typename EMBEDDING_SUBSYS_TYPE::CONST_STENCIL_PROXY_TYPE
                        (EMBEDDING_SUBSYS_TYPE::*)( int ) const,
                    &EMBEDDING_SUBSYS_TYPE::Stencil_Proxy
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
                As_Const_Array_View(ztaz_embedding_subsys.index_of_stencil_index),
                BOUND_FAST_MEM_FN<
                    typename ZTAZ_EMBEDDING_SUBSYS_TYPE::STENCIL_PROXY_TYPE
                        (ZTAZ_EMBEDDING_SUBSYS_TYPE::*)( int ),
                    &ZTAZ_EMBEDDING_SUBSYS_TYPE::Stencil_Proxy_Of_Stencil_Index
                >(ztaz_embedding_subsys)
            );
            std::cout << timer.Elapsed() << " s" << std::endl;

#if 0
            {
                std::cout << "Verifying correctness of Z^T*A*Z system via multiplication by standard basis vectors...";
                std::cout.flush();
                timer.Restart();
                ARRAY<T> e(n_index); // init'ed to 0
                ARRAY<T> ztaz_e_a(n_index, false); // uninit'ed
                ARRAY<T> ztaz_e_b(n_index, false); // uninit'ed
                T max_abs_diff = 0;
                for(int index = 1; index <= n_index; ++index) {
                    if(aggregate_constraint_system.indy_index_of_index.Contains(index))
                        continue;
                    e(index) = static_cast<T>(1);
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
                        const int index_indy = aggregate_constraint_system.index_of_indy_index(indy_index);
                        ztaz_e_b(index_indy) = static_cast<T>(0);
                    }
                    const T local_max_abs_diff = ARRAYS_COMPUTATIONS::Maxabs(ztaz_e_a - ztaz_e_b);
                    max_abs_diff = std::max(max_abs_diff, local_max_abs_diff);
                    e(index) = static_cast<T>(0);
                }
                std::cout << timer.Elapsed() << " s" << std::endl;
                std::cout << "  max abs diff = " << max_abs_diff << std::endl;
            }
#endif // #if 0|1

            if(main_params.general.randomized_check) {
                std::cout << "Verifying correctness of Z^T*A*Z system via multiplication by a random vector...";
                std::cout.flush();
                timer.Restart();
                ARRAY<T> u_random(n_index); // init'ed to 0
                for(int index = 1; index <= n_index; ++index)
                    if(
                        system.Stencil_N_Nonzero(index) != 0
                     && !aggregate_constraint_system.indy_index_of_index.Contains(index)
                    )
                        u_random(index) = rand();
                ARRAY<T> ztaz_u_random_a(n_index); // init'ed to 0
                {
                    ARRAY<T> z_u_random(u_random);
                    aggregate_constraint_system.Apply_Z(z_u_random);
                    system.Apply(z_u_random, ztaz_u_random_a);
                    aggregate_constraint_system.Apply_Z_Transpose(ztaz_u_random_a);
                }
                ARRAY<T> ztaz_u_random_b(n_index); // init'ed to 0
                regular_subsys.Apply(u_random, ztaz_u_random_b);
                ztaz_embedding_subsys.Apply(u_random, ztaz_u_random_b);
                for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
                    const int index = aggregate_constraint_system.index_of_indy_index(indy_index);
                    ztaz_u_random_b(index) = static_cast<T>(0);
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
                const int index = aggregate_constraint_system.index_of_indy_index(indy_index);
                const T c_value = aggregate_constraint_rhs(indy_index)
                                / aggregate_constraint_system.value_of_indy_index(indy_index);
                system.Apply_Transpose(index, ztaz_system_rhs, -c_value);
            }
            aggregate_constraint_system.Apply_Z_Transpose(ztaz_system_rhs);
            std::cout << timer.Elapsed() << " s" << std::endl;

            // Not necessary since all indys are virtual.
#if 0
            std::cout << "Zero'ing indy stencils in Z^T*A*Z regular subsystem...";
            std::cout.flush();
            timer.Restart();
            for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
                const int index = aggregate_constraint_system.index_of_indy_index(indy_index);
                assert(!ztaz_embedding_subsys.stencil_index_of_index.Contains(index));
                if(index_transform.Is_Material_Index(index)) {
                    const int linear_index = grid_index_of_index(index);
                    regular_subsys.Zero_Stencil(linear_index);
                }
            }
            std::cout << timer.Elapsed() << " s" << std::endl;
#endif // #if 0

            const ZTAZ_SYSTEM_TYPE ztaz_system(boost::fusion::make_vector(
                boost::ref(regular_subsys),
                boost::ref(ztaz_embedding_subsys))
            );

            std::cout << "Evaluating Z^T*A*Z residual norm of continuous solution...";
            std::cout.flush();
            timer.Restart();
            {
                ARRAY<T> ztaz_residual(-ztaz_system_rhs);
                ztaz_system.Apply(u_continuous, ztaz_residual);
                std::cout << timer.Elapsed() << " s" << std::endl;
                std::cout << "  " << ARRAYS_COMPUTATIONS::Maxabs(ztaz_residual) << std::endl;
            }

            switch(main_params.solver.solver_id) {
            case SOLVER_PARAMS::SOLVER_ID_NULL:
                return 0;
            case SOLVER_PARAMS::SOLVER_ID_PHYSBAM_CG:
                std::cout << "WARNING: Solver \"physbam-cg\" not yet implemented for Dirichlet problems." << std::endl;
                break;
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
                        main_params.general.n_thread,
                        ztaz_system, ztaz_system_rhs,
                        false, // has_constant_vectors_in_null_space
                        main_params.solver.max_iterations,
                        main_params.solver.relative_tolerance,
                        main_params.solver.absolute_tolerance,
                        main_params.solver.print_residuals,
                        main_params.solver.precondition,
                        u_approx
                    )
                ));
                std::cout << "[Solving with PETSc CG solver...] " << timer.Elapsed() << " s" << std::endl;
                break;
#endif // #ifdef PHYSBAM_NO_PETSC
            case SOLVER_PARAMS::SOLVER_ID_MG:
                std::cout << "WARNING: Solver \"mg\" not yet implemented for Neumann problems." << std::endl;
                break;
            case SOLVER_PARAMS::SOLVER_ID_MGPCG:
                std::cout << "WARNING: Solver \"mgpcg\" not yet implemented for Neumann problems." << std::endl;
                break;
            default:
                assert(false);
            }

            std::cout << "Evaluating c + Z*v...";
            std::cout.flush();
            timer.Restart();
            aggregate_constraint_system.Apply_Z(u_approx);
            for(int indy_index = 1; indy_index <= n_indy; ++indy_index) {
                const int index = aggregate_constraint_system.index_of_indy_index(indy_index);
                const T c_value = aggregate_constraint_rhs(indy_index)
                                / aggregate_constraint_system.value_of_indy_index(indy_index);
                u_approx(index) += c_value;
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
        problem.negative, main_params,
        -1,
        index_transform.sign_of_grid_index,
        Make_Compose_Function(
            As_Const_Array_View(regular_subsys.sign_of_cell_index),
            cell_multi_index_bound
        ),
        Make_Const_Array_View(multi_index_bound.Size(), u_approx.Get_Array_Pointer()),
        max_u_error, max_grad_u_error
    );
    Evaluate_Error(
        problem.positive, main_params,
        +1,
        index_transform.sign_of_grid_index,
        Make_Compose_Function(
            As_Const_Array_View(regular_subsys.sign_of_cell_index),
            cell_multi_index_bound
        ),
        Make_Const_Array_View(multi_index_bound.Size(), u_approx.Get_Array_Pointer()),
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
template int Build_And_Solve_Interface_System<T,D>( \
    const EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS& problem, \
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
