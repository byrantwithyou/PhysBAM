//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <iosfwd>

#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/mpl/vector/vector10.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/ref.hpp>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Functional/APPLY_MINUS_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/IF_ELSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Grid/ASSIGN_SIGN_TO_INDEX_GRID_VISITOR.h>
#include <Jeffrey_Utilities/Grid/Visit_Cells_With_Sign_Via_Fine_Vertex_Sign.h>
#include <Jeffrey_Utilities/Multi_Index/FINE_MULTI_INDEX_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_CUBE.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_X_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/Visit_Multi_Index_Box_Boundary.h>
#include <Jeffrey_Utilities/ONSTREAM.h>
#include <Jeffrey_Utilities/Stencils/INDEX_TRANSFORM_STENCIL_PROXY.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <Jeffrey_Utilities/VISITOR_SEQUENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "BETA_GRAD_U_DOT_N.h"
#include "Build_Embedding_Subsys.h"
#include "Build_Interface_Regular_Subsys.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "EMBEDDING_UNSTRUCTURED_SUBSYS.h"
#include "Init_Cell_Local_Embedding_Interface_System.h"
#include "INTERFACE_CONSTRAINT_SYSTEM.h"
#include "INTERFACE_INDEX_TRANSFORM.h"
#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "Print_System_Statistics.h"
#include "SET_DIRICHLET_GRID_BC_VISITOR.h"
#include "SYSTEM_SUM.h"

#include "Build_Interface_System.h"

namespace PhysBAM
{

template const int& HASHTABLE<int,int>::Get(const int&) const;

namespace Multigrid_Embedded_Poisson
{

namespace
{

template< class T, int D, class T_EMBEDDING_SUBSYS >
struct POST_EMBEDDING_INIT_VISITOR;
template< class T, int D, class T_EMBEDDING_SUBSYS >
struct EMBEDDING_CELL_VISITOR;

} // namespace

template< class T, int D, class T_EMBEDDING_SUBSYS >
int Build_Interface_System(
    const typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS& problem,
    const MAIN_PARAMS<T,D>& main_params,
    const ARRAY_VIEW<const T> phi_of_fine_index,
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys,
    T_EMBEDDING_SUBSYS& embedding_subsys,
    ARRAY<T>& system_rhs,
    INTERFACE_CONSTRAINT_SYSTEM<T,D>& constraint_system,
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
    assert(system_rhs.Size() == 0);

    typedef typename Result_Of::MAKE_COMPOSE_FUNCTION<
        SIGN_FUNCTION, ARRAY_VIEW<const T>, MULTI_INDEX_BOUND<D>
    >::type SIGN_OF_FINE_INDEX_TYPE;
    const SIGN_OF_FINE_INDEX_TYPE sign_of_fine_index =
        Make_Compose_Function(SIGN_FUNCTION(), phi_of_fine_index, fine_multi_index_bound);
    const MULTI_INDEX_X_FUNCTION< T, D, MULTI_INDEX_BOUND<D> >
        x_of_index(min_x, max_x, multi_index_bound);

    typedef SYSTEM_SUM< boost::mpl::vector2<
        DOMAIN_REGULAR_CROSS_SUBSYS<T,D>&,
        T_EMBEDDING_SUBSYS&
    > > SYSTEM_TYPE;
    SYSTEM_TYPE system(boost::fusion::make_vector(
        boost::ref(regular_subsys),
        boost::ref(embedding_subsys))
    );

    lout << "Initializing sign_of_cell_index...";
    lout.flush();
    timer.Restart();
    Visit_Cells_With_Sign_Via_Fine_Vertex_Sign_MT<2>(
        main_params.general.n_thread,
        cell_multi_index_bound,
        sign_of_fine_index,
        Make_Assign_Sign_To_Index_Grid_Visitor(
            Make_Array_Wrapper_Function(regular_subsys.sign_of_cell_index)
        ),
        -1 // sign_of_zero
    );
    lout << timer.Elapsed() << " s" << std::endl;

    {
        typedef BOUND_FAST_MEM_FN<
            const int& (HASHTABLE<int,int>::*)( const int& ) const,
            &HASHTABLE<int,int>::Get
        > VIRTUAL_INDEX_OFFSET_OF_GRID_INDEX_TYPE;
        typedef ARRAY_VIEW<const int> GRID_INDEX_OF_VIRTUAL_INDEX_OFFSET_TYPE;
        typedef typename Result_Of::MAKE_COMPOSE_FUNCTION<
            SIGN_OF_FINE_INDEX_TYPE,
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
                sign_of_fine_index,
                FINE_MULTI_INDEX_FUNCTION<2>(),
                multi_index_bound
            )
        );

        typedef BOUND_FAST_MEM_FN<
            int (INDEX_TRANSFORM_TYPE::*)( int ) const,
            &INDEX_TRANSFORM_TYPE::Grid_Index_Of_Index
        > GRID_INDEX_OF_INDEX_TYPE;
        typedef COMPOSE_FUNCTION<
            MULTI_INDEX_BOUND<D>,
            GRID_INDEX_OF_INDEX_TYPE
        > MULTI_INDEX_OF_INDEX_TYPE;
        const MULTI_INDEX_OF_INDEX_TYPE multi_index_of_index(
            multi_index_bound,
            GRID_INDEX_OF_INDEX_TYPE(index_transform)
        );

        typedef COMPOSE_FUNCTION<
            typename INDEX_TRANSFORM_TYPE::INDEX_OF_SIGNED_GRID_INDEX_FUNCTION,
            MULTI_INDEX_BOUND<D>
        > INDEX_OF_SIGNED_MULTI_INDEX_TYPE;
        const INDEX_OF_SIGNED_MULTI_INDEX_TYPE index_of_negative_multi_index(
            index_transform.Index_Of_Signed_Grid_Index_Function(-1),
            multi_index_bound
        );
        const INDEX_OF_SIGNED_MULTI_INDEX_TYPE index_of_positive_multi_index(
            index_transform.Index_Of_Signed_Grid_Index_Function(+1),
            multi_index_bound
        );

        typedef INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION<
            MULTI_INDEX_OF_INDEX_TYPE,
            INDEX_OF_SIGNED_MULTI_INDEX_TYPE
        > SIGNED_INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION_TYPE;
        const SIGNED_INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION_TYPE
            negative_index_transform_stencil_proxy_function(
                multi_index_of_index,
                index_of_negative_multi_index
            );
        const SIGNED_INDEX_TRANSFORM_STENCIL_PROXY_FUNCTION_TYPE
            positive_index_transform_stencil_proxy_function(
                multi_index_of_index,
                index_of_positive_multi_index
            );

        const BOUND_FAST_MEM_FN<
            typename T_EMBEDDING_SUBSYS::STENCIL_PROXY_TYPE (T_EMBEDDING_SUBSYS::*)( int ),
            &T_EMBEDDING_SUBSYS::Stencil_Proxy
        > system_stencil_proxy_of_index(embedding_subsys);

        typedef BOUND_FAST_MEM_FN<
            const int& (HASHTABLE<int,int>::*)( const int& ) const,
            &HASHTABLE<int,int>::Get
        > CONSTRAINT_STENCIL_INDEX_OF_CELL_LINEAR_INDEX_TYPE;
        typedef COMPOSE_FUNCTION<
            CONSTRAINT_STENCIL_INDEX_OF_CELL_LINEAR_INDEX_TYPE,
            MULTI_INDEX_BOUND<D>
        > CONSTRAINT_STENCIL_INDEX_OF_CELL_MULTI_INDEX_TYPE;
        const CONSTRAINT_STENCIL_INDEX_OF_CELL_MULTI_INDEX_TYPE
            constraint_stencil_index_of_cell_multi_index(
                CONSTRAINT_STENCIL_INDEX_OF_CELL_LINEAR_INDEX_TYPE(
                    constraint_system.stencil_index_of_cell_index
                ),
                cell_multi_index_bound
            );
        typedef BOUND_FAST_MEM_FN<
            typename INTERFACE_CONSTRAINT_SYSTEM<T,D>::STENCIL_PROXY_TYPE
                (INTERFACE_CONSTRAINT_SYSTEM<T,D>::*)( int ),
            &INTERFACE_CONSTRAINT_SYSTEM<T,D>::Stencil_Proxy
        > CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_STENCIL_INDEX_TYPE;
        const COMPOSE_FUNCTION<
            CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_STENCIL_INDEX_TYPE,
            CONSTRAINT_STENCIL_INDEX_OF_CELL_MULTI_INDEX_TYPE
        > constraint_stencil_proxy_of_cell_multi_index(
            CONSTRAINT_STENCIL_PROXY_OF_CONSTRAINT_STENCIL_INDEX_TYPE(constraint_system),
            constraint_stencil_index_of_cell_multi_index
        );

        Build_Embedding_Subsys(
            main_params.general.n_thread,
            multi_index_bound,
            As_Const_Array_View(regular_subsys.sign_of_cell_index),
            POST_EMBEDDING_INIT_VISITOR< T, D, T_EMBEDDING_SUBSYS >(
                multi_index_bound.Size(),
                embedding_subsys, system_rhs,
                constraint_system, constraint_rhs
            ),
            Make_Init_Cell_Local_Embedding_Interface_System_Visitor(
                min_x, max_x, multi_index_bound,
                Make_Compose_Function(phi_of_fine_index, fine_multi_index_bound),
                Make_Compose_Function(problem.negative.beta, x_of_index),
                Make_Compose_Function(problem.negative.f, x_of_index),
                Make_Compose_Function(problem.positive.beta, x_of_index),
                Make_Compose_Function(problem.positive.f, x_of_index),
                Make_Apply_Minus_Function(problem.positive.u, problem.negative.u),
                Make_Apply_Minus_Function(
                    Make_Beta_Grad_U_Dot_N(problem.positive.beta, problem.positive.grad_u),
                    Make_Beta_Grad_U_Dot_N(problem.negative.beta, problem.negative.grad_u)
                ),
                0.0f, // min_dist_to_vertex
                -1, // sign_of_zero
                Make_Compose_Function(
                    negative_index_transform_stencil_proxy_function,
                    system_stencil_proxy_of_index,
                    index_of_negative_multi_index
                ),
                Make_Compose_Function(
                    Make_Array_Wrapper_Function(system_rhs),
                    index_of_negative_multi_index
                ),
                Make_Compose_Function(
                    positive_index_transform_stencil_proxy_function,
                    system_stencil_proxy_of_index,
                    index_of_positive_multi_index
                ),
                Make_Compose_Function(
                    Make_Array_Wrapper_Function(system_rhs),
                    index_of_positive_multi_index
                ),
                Make_Compose_Function(
                    negative_index_transform_stencil_proxy_function,
                    constraint_stencil_proxy_of_cell_multi_index
                ),
                Make_Compose_Function(
                    positive_index_transform_stencil_proxy_function,
                    constraint_stencil_proxy_of_cell_multi_index
                ),
                Make_Compose_Function(
                    Make_Array_Wrapper_Function(constraint_rhs),
                    constraint_stencil_index_of_cell_multi_index
                )
            ),
            embedding_subsys.index_of_stencil_index,
            constraint_system.cell_index_of_stencil_index,
            lout
        );
        constraint_system.Init_Stencils_Containing_Index();
    }

    Build_Interface_Regular_Subsys(
        main_params.general.n_thread,
        dx, cell_multi_index_bound,
        Make_Compose_Function(problem.negative.beta, x_of_index),
        Make_Compose_Function(problem.negative.f, x_of_index),
        Make_Compose_Function(problem.positive.beta, x_of_index),
        Make_Compose_Function(problem.positive.f, x_of_index),
        regular_subsys, As_Array_View(system_rhs),
        lout
    );

    const int n_embedding = embedding_subsys.stencils.Size();
    const int n_index = multi_index_bound.Size() + n_embedding;

    Print_System_Statistics(system, n_index, lout);

    lout << "Setting Dirichlet grid bc's...";
    lout.flush();
    timer.Restart();
    Visit_Multi_Index_Box_Boundary(
        multi_index_bound,
        Make_If_Else_Function(
            Make_Compose_Function(
                Make_Equal_Function(-1),
                sign_of_fine_index,
                FINE_MULTI_INDEX_FUNCTION<2>()
            ),
            Make_Compose_Function(
                Make_Visitor_Sequence(
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        system,
                        Make_Compose_Function(problem.negative.u, x_of_index),
                        Make_Array_Wrapper_Function(system_rhs)
                    ),
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        constraint_system,
                        Make_Compose_Function(problem.negative.u, x_of_index),
                        Make_Array_Wrapper_Function(constraint_rhs)
                    )
                ),
                multi_index_bound
            ),
            Make_Compose_Function(
                Make_Visitor_Sequence(
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        system,
                        Make_Compose_Function(problem.positive.u, x_of_index),
                        Make_Array_Wrapper_Function(system_rhs)
                    ),
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        constraint_system,
                        Make_Compose_Function(problem.positive.u, x_of_index),
                        Make_Array_Wrapper_Function(constraint_rhs)
                    )
                ),
                multi_index_bound
            )
        )
    );
    lout << timer.Elapsed() << " s" << std::endl;

    return 0;
}

namespace
{

template< class T, int D, class T_EMBEDDING_SUBSYS >
struct POST_EMBEDDING_INIT_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        POST_EMBEDDING_INIT_VISITOR,
        (( /******/ int const, n_material ))
        (( typename T_EMBEDDING_SUBSYS&, embedding_subsys ))
        (( typename ARRAY<T>&, system_rhs ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( INTERFACE_CONSTRAINT_SYSTEM<T,D> )) &, constraint_system ))
        (( typename ARRAY<T>&, constraint_rhs ))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        const int n_virtual = embedding_subsys.index_of_stencil_index.Size();
        const int n_embedding = 2 * n_virtual;

        embedding_subsys.index_of_stencil_index.Preallocate(n_embedding);
        for(int i = 1; i <= n_virtual; ++i)
            embedding_subsys.index_of_stencil_index.Append(n_material + i);
        embedding_subsys.Init_Stencil_Index_Of_Index();
        embedding_subsys.stencils.Exact_Resize(n_embedding, false); // uninit'ed

        system_rhs.Exact_Resize(n_material + n_virtual); // init'ed to 0

        const int n_constraint = constraint_system.cell_index_of_stencil_index.Size();
        constraint_system.Init_Stencil_Index_Of_Cell_Index();
        constraint_system.stencils.Exact_Resize(n_constraint);
        constraint_system.stencils_containing_index.Initialize_New_Table(n_embedding);
        constraint_rhs.Exact_Resize(n_constraint);
    }
};

} // namespace

#define EXPLICIT_INSTANTIATION( T, D ) \
    EXPLICIT_INSTANTIATION_HELPER( T, D, ( EMBEDDING_UNSTRUCTURED_SUBSYS<T> ) )
#define EXPLICIT_INSTANTIATION_HELPER( T, D, T_EMBEDDING_SUBSYS ) \
template int \
Build_Interface_System< T, D, BOOST_PP_SEQ_ENUM( T_EMBEDDING_SUBSYS ) >( \
    EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS const & problem, \
    const MAIN_PARAMS<T,D>& main_params, \
    const ARRAY_VIEW<const T> phi_of_fine_index, \
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys, \
    BOOST_PP_SEQ_ENUM( T_EMBEDDING_SUBSYS )& embedding_subsys, \
    ARRAY<T>& system_rhs, \
    INTERFACE_CONSTRAINT_SYSTEM<T,D>& constraint_system, \
    ARRAY<T>& constraint_rhs, \
    std::ostream& lout);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION
#undef EXPLICIT_INSTANTIATION_HELPER

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
