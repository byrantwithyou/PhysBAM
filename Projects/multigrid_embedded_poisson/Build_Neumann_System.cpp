//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <cassert>

#include <iostream>
#include <limits>

#include <boost/preprocessor/seq/enum.hpp>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/DIRECT_INIT_CTOR.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/BOUND_FAST_MEM_FN.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/SIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/VISIT_IF.h>
#include <Jeffrey_Utilities/Grid/ASSIGN_SIGN_TO_INDEX_GRID_VISITOR.h>
#include <Jeffrey_Utilities/Grid/Visit_Cells_With_Sign_Via_Fine_Vertex_Sign.h>
#include <Jeffrey_Utilities/IDENTITY_TYPE.h>
#include <Jeffrey_Utilities/Multi_Index/FINE_MULTI_INDEX_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_X_FUNCTION.h>
#include <Jeffrey_Utilities/Multi_Index/Visit_Multi_Index_Box_Boundary.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "BETA_GRAD_U_DOT_N.h"
#include "Build_Domain_Regular_Subsys.h"
#include "Build_Embedding_Subsys.h"
#include "DOMAIN_EMBEDDING_CUBE_SUBSYS.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "DOMAIN_SYSTEM.h"
#include "Init_Cell_Local_Embedding_Neumann_System.h"
#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "SET_DIRICHLET_GRID_BC_VISITOR.h"

#include "Build_Neumann_System.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace
{

template< class T_EMBEDDING_SUBSYS >
struct POST_EMBEDDING_INIT_VISITOR;
template< class T, int D, class T_EMBEDDING_SUBSYS >
struct EMBEDDING_CELL_VISITOR;

} // namespace

template< class T, int D, class T_EMBEDDING_SUBSYS >
int Build_Neumann_System(
    typename EXAMPLE_PARAMS<T,D>::NEUMANN_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    const ARRAY_VIEW<const T> phi_of_fine_index,
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys,
    T_EMBEDDING_SUBSYS& embedding_subsys,
    ARRAY_VIEW<T> system_rhs)
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

    DOMAIN_SYSTEM<
        DOMAIN_REGULAR_CROSS_SUBSYS<T,D>&,
        T_EMBEDDING_SUBSYS&
    > system(regular_subsys, embedding_subsys);

    std::cout << "Initializing sign_of_cell_index...";
    std::cout.flush();
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
        POST_EMBEDDING_INIT_VISITOR< T_EMBEDDING_SUBSYS >(embedding_subsys),
        EMBEDDING_CELL_VISITOR< T, D, T_EMBEDDING_SUBSYS >(
            problem, main_params,
            phi_of_fine_index,
            embedding_subsys, system_rhs
        ),
        embedding_subsys.linear_index_of_stencil_index,
        std::cout
    );

    Build_Domain_Regular_Subsys(
        main_params.general.n_thread,
        dx, cell_multi_index_bound,
        Make_Compose_Function(
            problem.beta,
            Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
        ),
        Make_Compose_Function(
            problem.f,
            Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
        ),
        regular_subsys, system_rhs,
        std::cout
    );

    {
        std::cout << "Computing system statistics...";
        std::cout.flush();
        timer.Restart();
        T max_diag = 0;
        T min_diag = std::numeric_limits<T>::infinity();
        T max_abs_stencil_sum = 0;
        for(int linear_index = 1; linear_index <= multi_index_bound.Size(); ++linear_index) {
            const T diag = system.Diag(linear_index);
            const T stencil_sum = system.Stencil_Sum(linear_index);
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
            Make_Set_Dirichlet_Grid_BC_Visitor(
                system,
                Make_Compose_Function(
                    problem.u,
                    Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
                ),
                Make_Array_Wrapper_Function(system_rhs)
            )
        )
    );
    std::cout << timer.Elapsed() << " s" << std::endl;

    return 0;
}

namespace
{

template< class T_EMBEDDING_SUBSYS >
struct POST_EMBEDDING_INIT_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        POST_EMBEDDING_INIT_VISITOR,
        (( typename T_EMBEDDING_SUBSYS&, embedding_subsys ))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        const int n_embedding = embedding_subsys.linear_index_of_stencil_index.Size();
        embedding_subsys.Init_Stencil_Index_Of_Linear_Index();
        embedding_subsys.stencils.Exact_Resize(n_embedding, false); // uninit'ed
        embedding_subsys.Zero_Stencils();
    }
};

template< class T, int D, class T_EMBEDDING_SUBSYS >
struct EMBEDDING_CELL_VISITOR
{
    typedef typename EXAMPLE_PARAMS<T,D>::NEUMANN_PARAMS NEUMANN_PARAMS_TYPE;
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        EMBEDDING_CELL_VISITOR,
        (( typename NEUMANN_PARAMS_TYPE const &, problem ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( MAIN_PARAMS<T,D> )) const &, main_params ))
        (( typename ARRAY_VIEW<const T> const, phi_of_fine_index ))
        (( typename T_EMBEDDING_SUBSYS&, embedding_subsys ))
        (( typename ARRAY_VIEW<T>&, system_rhs ))
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
        Init_Cell_Local_Embedding_Neumann_System(
            min_x, max_x, multi_index_bound,
            -1, // domain_sign
            Make_Compose_Function(phi_of_fine_index, fine_multi_index_bound),
            Make_Compose_Function(
                problem.beta,
                Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
            ),
            Make_Compose_Function(
                problem.f,
                Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
            ),
            Make_Beta_Grad_U_Dot_N(problem.beta, problem.grad_u),
            0.0f, // min_dist_to_vertex
            -1, // sign_of_zero
            cell_multi_index,
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
            )
        );
    }
};

} // namespace

#define EXPLICIT_INSTANTIATION( T, D ) \
    EXPLICIT_INSTANTIATION_HELPER( T, D, ( DOMAIN_EMBEDDING_CUBE_SUBSYS< T ) ( D > ) )
#define EXPLICIT_INSTANTIATION_HELPER( T, D, T_EMBEDDING_SUBSYS ) \
template int \
Build_Neumann_System< T, D, BOOST_PP_SEQ_ENUM( T_EMBEDDING_SUBSYS ) >( \
    EXAMPLE_PARAMS<T,D>::NEUMANN_PARAMS const & problem, \
    const MAIN_PARAMS<T,D>& main_params, \
    const ARRAY_VIEW<const T> phi_of_fine_index, \
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys, \
    BOOST_PP_SEQ_ENUM( T_EMBEDDING_SUBSYS )& embedding_subsys, \
    ARRAY_VIEW<T> system_rhs);
EXPLICIT_INSTANTIATION( float, 2 )
EXPLICIT_INSTANTIATION( float, 3 )
EXPLICIT_INSTANTIATION( double, 2 )
EXPLICIT_INSTANTIATION( double, 3 )
#undef EXPLICIT_INSTANTIATION
#undef EXPLICIT_INSTANTIATION_HELPER

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM
