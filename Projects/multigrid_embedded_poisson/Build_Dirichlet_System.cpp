//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <iosfwd>

#include <boost/preprocessor/seq/enum.hpp>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Functional/ARGUMENT_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/CONSTANT_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/STATIC_CAST_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/VISIT_IF.h>
#include <Jeffrey_Utilities/Grid/Cell_Sign_Via_Fine_Vertex_Sign.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Multi_Index/FINE_MULTI_INDEX_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_X_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/Visit_Multi_Index_Box_Boundary.h>
#include <Jeffrey_Utilities/ONSTREAM.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <Jeffrey_Utilities/VISITOR_SEQUENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "BETA_GRAD_U_DOT_N.h"
#include "Build_Domain_Regular_Subsys.h"
#include "Build_Embedding_Subsys.h"
#include "DIRICHLET_CONSTRAINT_SYSTEM.h"
#include "DOMAIN_EMBEDDING_CUBE_SUBSYS.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "DOMAIN_SYSTEM.h"
#include "Init_Cell_Local_Embedding_Dirichlet_System.h"
#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "POST_EMBEDDING_INIT_DIRICHLET_CONSTRAINT_VISITOR.h"
#include "POST_EMBEDDING_INIT_DOMAIN_VISITOR.h"
#include "Print_System_Statistics.h"
#include "SET_DIRICHLET_GRID_BC_VISITOR.h"
#include "SET_NEUMANN_OFFSET_GRID_BC_VISITOR.h"

#include "Build_Dirichlet_System.h"

namespace PhysBAM
{

template const int& HASHTABLE<int,int>::Get(const int&) const;

namespace Multigrid_Embedded_Poisson
{

template< class T, int D, class T_SIGN, class T_EMBEDDING_SUBSYS >
int Build_Dirichlet_System(
    const typename EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS& problem,
    const MAIN_PARAMS<T,D>& main_params,
    const ARRAY_VIEW<const T> phi_of_fine_index,
    const ARRAY_VIEW<const T_SIGN> sign_of_cell_index,
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys,
    T_EMBEDDING_SUBSYS& embedding_subsys,
    ARRAY_VIEW<T> system_rhs,
    DIRICHLET_CONSTRAINT_SYSTEM<T,D>& constraint_system,
    ARRAY<T>& constraint_rhs,
    std::ostream& lout = PhysBAM::nout)
{
    typedef VECTOR<int,D> MULTI_INDEX_TYPE;

    BASIC_TIMER timer;

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
    const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;

    const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
    const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);
    const VECTOR<T,D> dx = (max_x - min_x) / cell_multi_index_bound.max_multi_index;

    assert(phi_of_fine_index.Size() == fine_multi_index_bound.Size());
    assert(system_rhs.Size() == multi_index_bound.Size());

    typename Result_Of::MAKE_COMPOSE_FUNCTION<
        SIGN_FUNCTION, ARRAY_VIEW<const T>, MULTI_INDEX_BOUND<D>
    >::type const sign_of_fine_index = Make_Compose_Function(
        SIGN_FUNCTION(), phi_of_fine_index, fine_multi_index_bound
    );
    const MULTI_INDEX_X_FUNCTION< T, D, MULTI_INDEX_BOUND<D> >
        x_of_index(min_x, max_x, multi_index_bound);

#ifndef NDEBUG
    for(int cell_linear_index = 1; cell_linear_index <= cell_multi_index_bound.Size(); ++cell_linear_index) {
        const MULTI_INDEX_TYPE cell_multi_index = cell_multi_index_bound.Multi_Index(cell_linear_index);
        const int cell_sign1 = Cell_Sign_Via_Fine_Vertex_Sign<2>(sign_of_fine_index, cell_multi_index, -1);
        const int cell_sign2 = sign_of_cell_index(cell_linear_index);
        assert(cell_sign1 == cell_sign2);
    }
#endif // #ifndef NDEBUG

    DOMAIN_SYSTEM<
        DOMAIN_REGULAR_CROSS_SUBSYS<T,D>&,
        T_EMBEDDING_SUBSYS&
    > system(regular_subsys, embedding_subsys);

    {
        typedef BOUND_FAST_MEM_FN<
            const int& (HASHTABLE<int,int>::*)( const int& ) const,
            &HASHTABLE<int,int>::Get
        > CONSTRAINT_STENCIL_INDEX_OF_CELL_LINEAR_INDEX_TYPE;
        const COMPOSE_FUNCTION<
            CONSTRAINT_STENCIL_INDEX_OF_CELL_LINEAR_INDEX_TYPE,
            MULTI_INDEX_BOUND<D>
        > constraint_stencil_index_of_cell_multi_index(
            CONSTRAINT_STENCIL_INDEX_OF_CELL_LINEAR_INDEX_TYPE(
                constraint_system.stencil_index_of_cell_linear_index
            ),
            cell_multi_index_bound
        );

        Build_Embedding_Subsys(
            main_params.general.n_thread,
            multi_index_bound,
            sign_of_cell_index,
            Make_Visitor_Sequence(
                Make_Post_Embedding_Init_Domain_Visitor(embedding_subsys),
                POST_EMBEDDING_INIT_DIRICHLET_CONSTRAINT_VISITOR<T,D>(
                    constraint_system, constraint_rhs
                )
            ),
            Make_Init_Cell_Local_Embedding_Dirichlet_System_Visitor(
                min_x, max_x, multi_index_bound,
                -1, // domain_sign
                Make_Compose_Function(phi_of_fine_index, fine_multi_index_bound),
                Make_Compose_Function(problem.beta, x_of_index),
                Make_Compose_Function(problem.f, x_of_index),
                Make_Constant_Function(problem.u),
                0.0f, // min_dist_to_vertex
                -1, // sign_of_zero
                Make_Compose_Function(
                    BOUND_FAST_MEM_FN<
                        typename T_EMBEDDING_SUBSYS::MULTI_INDEX_STENCIL_PROXY_TYPE
                            (T_EMBEDDING_SUBSYS::*)( int ),
                        &T_EMBEDDING_SUBSYS::Multi_Index_Stencil_Proxy
                    >(embedding_subsys),
                    multi_index_bound
                ),
                Make_Compose_Function(
                    Make_Array_Wrapper_Function(system_rhs),
                    multi_index_bound
                ),
                Make_Compose_Function(
                    BOUND_FAST_MEM_FN<
                        typename DIRICHLET_CONSTRAINT_SYSTEM<T,D>::MULTI_INDEX_STENCIL_PROXY_TYPE
                            (DIRICHLET_CONSTRAINT_SYSTEM<T,D>::*)( int ),
                        &DIRICHLET_CONSTRAINT_SYSTEM<T,D>::Multi_Index_Stencil_Proxy
                    >(constraint_system),
                    constraint_stencil_index_of_cell_multi_index
                ),
                Make_Compose_Function(
                    Make_Array_Wrapper_Function(constraint_rhs),
                    constraint_stencil_index_of_cell_multi_index
                )
            ),
            embedding_subsys.linear_index_of_stencil_index,
            constraint_system.cell_linear_index_of_stencil_index,
            lout
        );
    }

    Build_Domain_Regular_Subsys(
        main_params.general.n_thread,
        dx, cell_multi_index_bound,
        sign_of_cell_index,
        Make_Compose_Function(problem.beta, x_of_index),
        Make_Compose_Function(problem.f, x_of_index),
        regular_subsys, system_rhs,
        lout
    );

    Print_System_Statistics(system, multi_index_bound.Size(), lout);

    switch(main_params.example.grid_bc_id) {
    case EXAMPLE_PARAMS_BASE::GRID_BC_ID_NEUMANN_OFFSET:
#if 0
        lout << "Setting Neumann offset grid bc's...";
        lout.flush();
        timer.Restart();
        Visit_Multi_Index_Box_Boundary(
            multi_index_bound,
            Make_Visit_If(
                Make_Compose_Function(
                    Make_Equal_Function(-1),
                    sign_of_fine_index,
                    FINE_MULTI_INDEX_FUNCTION<2>()
                ),
                Make_Set_Neumann_Offset_Grid_BC_Visitor(
                    dx, multi_index_bound,
                    Make_Compose2_Function(
                        Make_Beta_Grad_U_Dot_N(problem.beta, problem.grad_u),
                        Make_Compose_Function(
                            Make_Multi_Index_X_Function(min_x, max_x, fine_multi_index_bound),
                            ARGUMENT_FUNCTION<1>()
                        ),
                        Make_Compose_Function(
                            STATIC_CAST_FUNCTION< VECTOR<T,D> >(),
                            ARGUMENT_FUNCTION<2>()
                        )
                    ),
                    Make_Compose_Function(
                        Make_Array_Wrapper_Function(system_rhs),
                        multi_index_bound
                    )
                )
            )
        );
        lout << timer.Elapsed() << " s" << std::endl;
        break;
#endif
    case EXAMPLE_PARAMS_BASE::GRID_BC_ID_DIRICHLET:
        lout << "Setting Dirichlet grid bc's...";
        lout.flush();
        timer.Restart();
        Visit_Multi_Index_Box_Boundary(
            multi_index_bound,
            Make_Visit_If(
                Make_Compose_Function(
                    Make_Equal_Function(-1),
                    SIGN_FUNCTION(),
                    phi_of_fine_index,
                    fine_multi_index_bound,
                    FINE_MULTI_INDEX_FUNCTION<2>()
                ),
                Make_Visitor_Sequence(
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        system,
                        Make_Compose_Function(problem.u, x_of_index),
                        Make_Array_Wrapper_Function(system_rhs)
                    ),
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        constraint_system,
                        Make_Compose_Function(problem.u, x_of_index),
                        Make_Array_Wrapper_Function(constraint_rhs)
                    )
                )
            )
        );
        lout << timer.Elapsed() << " s" << std::endl;
        break;
    default:
        assert(false);
    }

    return 0;
}

#define EXPLICIT_INSTANTIATION( T, D, T_SIGN ) \
    EXPLICIT_INSTANTIATION_HELPER( T, D, T_SIGN, ( DOMAIN_EMBEDDING_CUBE_SUBSYS< T ) ( D > ) )
#define EXPLICIT_INSTANTIATION_HELPER( T, D, T_SIGN, T_EMBEDDING_SUBSYS ) \
template int \
Build_Dirichlet_System< T, D, T_SIGN, BOOST_PP_SEQ_ENUM( T_EMBEDDING_SUBSYS ) >( \
    EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS const & problem, \
    const MAIN_PARAMS<T,D>& main_params, \
    const ARRAY_VIEW<const T> phi_of_fine_index, \
    const ARRAY_VIEW<const T_SIGN> sign_of_cell_index, \
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys, \
    BOOST_PP_SEQ_ENUM( T_EMBEDDING_SUBSYS )& embedding_subsys, \
    ARRAY_VIEW<T> system_rhs, \
    DIRICHLET_CONSTRAINT_SYSTEM<T,D>& constraint_system, \
    ARRAY<T>& constraint_rhs, \
    std::ostream& lout);
EXPLICIT_INSTANTIATION( float, 2, signed char )
EXPLICIT_INSTANTIATION( float, 3, signed char )
EXPLICIT_INSTANTIATION( double, 2, signed char )
EXPLICIT_INSTANTIATION( double, 3, signed char )
#undef EXPLICIT_INSTANTIATION
#undef EXPLICIT_INSTANTIATION_HELPER

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
