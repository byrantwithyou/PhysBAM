//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <iostream>
#include <limits>

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
    ARRAY<T>& constraint_rhs)
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

    typedef SYSTEM_SUM< boost::mpl::vector2<
        DOMAIN_REGULAR_CROSS_SUBSYS<T,D>&,
        T_EMBEDDING_SUBSYS&
    > > SYSTEM_TYPE;
    SYSTEM_TYPE system(boost::fusion::make_vector(
        boost::ref(regular_subsys),
        boost::ref(embedding_subsys))
    );

    std::cout << "Initializing sign_of_cell_index...";
    std::cout.flush();
    timer.Restart();
    Visit_Cells_With_Sign_Via_Fine_Vertex_Sign_MT<2>(
        main_params.general.n_thread,
        cell_multi_index_bound,
        Make_Compose_Function(
            SIGN_FUNCTION(),
            phi_of_fine_index,
            fine_multi_index_bound
        ),
        Make_Assign_Sign_To_Index_Grid_Visitor(
            Make_Array_Wrapper_Function(regular_subsys.sign_of_cell_index)
        ),
        -1 // sign_of_zero
    );
    std::cout << timer.Elapsed() << " s" << std::endl;

    Build_Embedding_Subsys(
        main_params.general.n_thread,
        multi_index_bound,
        As_Const_Array_View(regular_subsys.sign_of_cell_index),
        POST_EMBEDDING_INIT_VISITOR< T, D, T_EMBEDDING_SUBSYS >(
            main_params,
            embedding_subsys, system_rhs,
            constraint_system, constraint_rhs
        ),
        EMBEDDING_CELL_VISITOR< T, D, T_EMBEDDING_SUBSYS >(
            problem, main_params,
            phi_of_fine_index,
            embedding_subsys, system_rhs,
            constraint_system, constraint_rhs
        ),
        embedding_subsys.index_of_stencil_index,
        constraint_system.cell_index_of_stencil_index,
        std::cout
    );
    constraint_system.Init_Stencils_Containing_Index();

    Build_Interface_Regular_Subsys(
        main_params.general.n_thread,
        dx, cell_multi_index_bound,
        Make_Compose_Function(
            problem.negative.beta,
            Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
        ),
        Make_Compose_Function(
            problem.positive.beta,
            Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
        ),
        Make_Compose_Function(
            problem.negative.f,
            Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
        ),
        Make_Compose_Function(
            problem.positive.f,
            Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
        ),
        regular_subsys, As_Array_View(system_rhs),
        std::cout
    );

    const int n_embedding = embedding_subsys.stencils.Size();
    const int n_index = multi_index_bound.Size() + n_embedding;

    {
        std::cout << "Computing system statistics...";
        std::cout.flush();
        timer.Restart();
        T max_diag = 0;
        T min_diag = std::numeric_limits<T>::infinity();
        T max_abs_stencil_sum = 0;
        for(int index = 1; index <= n_index; ++index) {
            const T diag = system.Diag(index);
            const T stencil_sum = system.Stencil_Sum(index);
            if(diag != 0) {
                max_diag = std::max(max_diag, diag);
                min_diag = std::min(min_diag, diag);
            }
            max_abs_stencil_sum = std::max(max_abs_stencil_sum, std::abs(stencil_sum));
        }
        std::cout << timer.Elapsed() << " s" << std::endl;
        std::cout << "  max diag = " << max_diag << '\n'
                  << "  min diag = " << min_diag << '\n'
                  << "    ratio = " << max_diag / min_diag << '\n'
                  << "  max abs stencil sum = " << max_abs_stencil_sum
                  << std::endl;
    }

    std::cout << "Setting Dirichlet grid bc's...";
    std::cout.flush();
    timer.Restart();
    typename Result_Of::MAKE_COMPOSE_FUNCTION<
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE,
        MULTI_INDEX_X_FUNCTION< T, D, MULTI_INDEX_BOUND<D> >
    >::type u_negative_of_index = Make_Compose_Function(
        problem.negative.u,
        Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
    );
    typename Result_Of::MAKE_COMPOSE_FUNCTION<
        typename EXAMPLE_PARAMS<T,D>::SCALAR_FUNCTION_TYPE,
        MULTI_INDEX_X_FUNCTION< T, D, MULTI_INDEX_BOUND<D> >
    >::type u_positive_of_index = Make_Compose_Function(
        problem.positive.u,
        Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
    );
    Visit_Multi_Index_Box_Boundary(
        multi_index_bound,
        Make_If_Else_Function(
            Make_Compose_Function(
                Make_Equal_Function(-1),
                SIGN_FUNCTION(),
                phi_of_fine_index,
                fine_multi_index_bound,
                FINE_MULTI_INDEX_FUNCTION<2>()
            ),
            Make_Compose_Function(
                Make_Visitor_Sequence(
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        system,
                        u_negative_of_index,
                        Make_Array_Wrapper_Function(system_rhs)
                    ),
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        constraint_system,
                        u_negative_of_index,
                        Make_Array_Wrapper_Function(constraint_rhs)
                    )
                ),
                multi_index_bound
            ),
            Make_Compose_Function(
                Make_Visitor_Sequence(
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        system,
                        u_positive_of_index,
                        Make_Array_Wrapper_Function(system_rhs)
                    ),
                    Make_Set_Dirichlet_Grid_BC_Visitor(
                        constraint_system,
                        u_positive_of_index,
                        Make_Array_Wrapper_Function(constraint_rhs)
                    )
                ),
                multi_index_bound
            )
        )
    );
    std::cout << timer.Elapsed() << " s" << std::endl;

    return 0;
}

namespace
{

template< class T, int D, class T_EMBEDDING_SUBSYS >
struct POST_EMBEDDING_INIT_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        POST_EMBEDDING_INIT_VISITOR,
        (( typename typename PHYSBAM_IDENTITY_TYPE(( MAIN_PARAMS<T,D> )) const &, main_params ))
        (( typename T_EMBEDDING_SUBSYS&, embedding_subsys ))
        (( typename ARRAY<T>&, system_rhs ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( INTERFACE_CONSTRAINT_SYSTEM<T,D> )) &, constraint_system ))
        (( typename ARRAY<T>&, constraint_rhs ))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
        const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;

        const int n_virtual = embedding_subsys.index_of_stencil_index.Size();
        const int n_embedding = 2 * n_virtual;
        embedding_subsys.index_of_stencil_index.Preallocate(n_embedding);
        for(int i = 1; i <= n_virtual; ++i)
            embedding_subsys.index_of_stencil_index.Append(multi_index_bound.Size() + i);
        embedding_subsys.Init_Stencil_Index_Of_Index();
        embedding_subsys.stencils.Exact_Resize(n_embedding, false); // uninit'ed

        system_rhs.Exact_Resize(multi_index_bound.Size() + n_virtual); // init'ed to 0

        const int n_constraint = constraint_system.cell_index_of_stencil_index.Size();
        constraint_system.Init_Stencil_Index_Of_Cell_Index();
        constraint_system.stencils.Exact_Resize(n_constraint);
        constraint_system.stencils_containing_index.Initialize_New_Table(n_embedding);
        constraint_rhs.Exact_Resize(n_constraint);
    }
};

template< class T, int D, class T_EMBEDDING_SUBSYS >
struct EMBEDDING_CELL_VISITOR
{
    typedef typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS INTERFACE_PARAMS_TYPE;
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        EMBEDDING_CELL_VISITOR,
        (( typename INTERFACE_PARAMS_TYPE const &, problem ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( MAIN_PARAMS<T,D> )) const &, main_params ))
        (( typename ARRAY_VIEW<const T> const, phi_of_fine_index ))
        (( typename T_EMBEDDING_SUBSYS&, embedding_subsys ))
        (( typename ARRAY<T>&, system_rhs ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( INTERFACE_CONSTRAINT_SYSTEM<T,D> )) &, constraint_system ))
        (( typename ARRAY<T>&, constraint_rhs ))
    )

public:

    typedef void result_type;

    void operator()(const int cell_linear_index) const
    {
        typedef VECTOR<int,D> MULTI_INDEX_TYPE;

        const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
        const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
        const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;
        const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
        const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);
        const MULTI_INDEX_TYPE cell_multi_index = cell_multi_index_bound.Multi_Index(cell_linear_index);
        const int constraint_stencil_index = constraint_system.stencil_index_of_cell_index.Get(cell_linear_index);

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

        typename Result_Of::MAKE_COMPOSE_FUNCTION<
            MULTI_INDEX_BOUND<D>,
            BOUND_FAST_MEM_FN<
                int (INDEX_TRANSFORM_TYPE::*)( int ) const,
                &INDEX_TRANSFORM_TYPE::Grid_Index_Of_Index
            >
        >::type const multi_index_of_index = Make_Compose_Function(
            multi_index_bound,
            PHYSBAM_BOUND_FAST_MEM_FN_TEMPLATE(
                index_transform,
                &INDEX_TRANSFORM_TYPE::Grid_Index_Of_Index
            )
        );

        typename Result_Of::MAKE_COMPOSE_FUNCTION<
            typename INDEX_TRANSFORM_TYPE::INDEX_OF_SIGNED_GRID_INDEX_FUNCTION,
            MULTI_INDEX_BOUND<D>
        >::type const index_of_negative_multi_index = Make_Compose_Function(
            index_transform.Index_Of_Signed_Grid_Index_Function(-1),
            multi_index_bound
        );
        typename Result_Of::MAKE_COMPOSE_FUNCTION<
            typename INDEX_TRANSFORM_TYPE::INDEX_OF_SIGNED_GRID_INDEX_FUNCTION,
            MULTI_INDEX_BOUND<D>
        >::type const index_of_positive_multi_index = Make_Compose_Function(
            index_transform.Index_Of_Signed_Grid_Index_Function(+1),
            multi_index_bound
        );

        Init_Cell_Local_Embedding_Interface_System(
            min_x, max_x, multi_index_bound,
            Make_Compose_Function(phi_of_fine_index, fine_multi_index_bound),
            Make_Compose_Function(
                problem.negative.beta,
                Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
            ),
            Make_Compose_Function(
                problem.negative.f,
                Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
            ),
            Make_Compose_Function(
                problem.positive.beta,
                Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
            ),
            Make_Compose_Function(
                problem.positive.f,
                Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
            ),
            Make_Apply_Minus_Function(problem.positive.u, problem.negative.u),
            Make_Apply_Minus_Function(
                Make_Beta_Grad_U_Dot_N(problem.positive.beta, problem.positive.grad_u),
                Make_Beta_Grad_U_Dot_N(problem.negative.beta, problem.negative.grad_u)
            ),
            0.0f, // min_dist_to_vertex
            -1, // sign_of_zero
            cell_multi_index,
            Make_Compose_Function(
                Make_Index_Transform_Stencil_Proxy_Function(
                    multi_index_of_index,
                    index_of_negative_multi_index
                ),
                BOUND_FAST_MEM_FN<
                    typename T_EMBEDDING_SUBSYS::STENCIL_PROXY_TYPE
                        (T_EMBEDDING_SUBSYS::*)( int ),
                    &T_EMBEDDING_SUBSYS::Stencil_Proxy
                >(embedding_subsys),
                index_of_negative_multi_index
            ),
            Make_Compose_Function(
                Make_Index_Transform_Stencil_Proxy_Function(
                    multi_index_of_index,
                    index_of_positive_multi_index
                ),
                BOUND_FAST_MEM_FN<
                    typename T_EMBEDDING_SUBSYS::STENCIL_PROXY_TYPE
                        (T_EMBEDDING_SUBSYS::*)( int ),
                    &T_EMBEDDING_SUBSYS::Stencil_Proxy
                >(embedding_subsys),
                index_of_positive_multi_index
            ),
            Make_Compose_Function(
                Make_Array_Wrapper_Function(system_rhs),
                index_of_negative_multi_index
            ),
            Make_Compose_Function(
                Make_Array_Wrapper_Function(system_rhs),
                index_of_positive_multi_index
            ),
            Make_Index_Transform_Stencil_Proxy(
                constraint_system.Stencil_Proxy(constraint_stencil_index),
                multi_index_of_index,
                index_of_negative_multi_index
            ),
            Make_Index_Transform_Stencil_Proxy(
                constraint_system.Stencil_Proxy(constraint_stencil_index),
                multi_index_of_index,
                index_of_positive_multi_index
            ),
            constraint_rhs(constraint_stencil_index)
        );
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
    ARRAY<T>& constraint_rhs);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION
#undef EXPLICIT_INSTANTIATION_HELPER

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
