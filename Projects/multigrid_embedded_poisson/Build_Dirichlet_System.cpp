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
#include <Jeffrey_Utilities/ONSTREAM.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <Jeffrey_Utilities/VISITOR_SEQUENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Build_Domain_Regular_Subsys.h"
#include "Build_Embedding_Subsys.h"
#include "DIRICHLET_CONSTRAINT_SYSTEM.h"
#include "DOMAIN_EMBEDDING_CUBE_SUBSYS.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "DOMAIN_SYSTEM.h"
#include "Init_Cell_Local_Embedding_Dirichlet_System.h"
#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "Print_System_Statistics.h"
#include "SET_DIRICHLET_GRID_BC_VISITOR.h"

#include "Build_Dirichlet_System.h"

namespace PhysBAM
{

template const int& HASHTABLE<int,int>::Get(const int&) const;

namespace Multigrid_Embedded_Poisson
{

namespace
{

template< class T, int D, class T_EMBEDDING_SUBSYS >
struct POST_EMBEDDING_INIT_VISITOR;

} // namespace

template< class T, int D, class T_EMBEDDING_SUBSYS >
int Build_Dirichlet_System(
    const typename EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS& problem,
    const MAIN_PARAMS<T,D>& main_params,
    const ARRAY_VIEW<const T> phi_of_fine_index,
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

    DOMAIN_SYSTEM<
        DOMAIN_REGULAR_CROSS_SUBSYS<T,D>&,
        T_EMBEDDING_SUBSYS&
    > system(regular_subsys, embedding_subsys);

    lout << "Initializing sign_of_cell_index...";
    lout.flush();
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
            As_Const_Array_View(regular_subsys.sign_of_cell_index),
            POST_EMBEDDING_INIT_VISITOR< T, D, T_EMBEDDING_SUBSYS >(
                embedding_subsys,
                constraint_system, constraint_rhs
            ),
            Make_Init_Cell_Local_Embedding_Dirichlet_System_Visitor(
                min_x, max_x, multi_index_bound,
                -1, // domain_sign
                Make_Compose_Function(phi_of_fine_index, fine_multi_index_bound),
                Make_Compose_Function(problem.beta, x_of_index),
                Make_Compose_Function(problem.f, x_of_index),
                problem.u,
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
        Make_Compose_Function(problem.beta, x_of_index),
        Make_Compose_Function(problem.f, x_of_index),
        regular_subsys, system_rhs,
        lout
    );

    Print_System_Statistics(system, multi_index_bound.Size(), lout);

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

    return 0;
}

namespace
{

template< class T, int D, class T_EMBEDDING_SUBSYS >
struct POST_EMBEDDING_INIT_VISITOR
{
    PHYSBAM_DIRECT_INIT_CTOR_DECLARE_PRIVATE_MEMBERS(
        POST_EMBEDDING_INIT_VISITOR,
        (( typename T_EMBEDDING_SUBSYS&, embedding_subsys ))
        (( typename typename PHYSBAM_IDENTITY_TYPE(( DIRICHLET_CONSTRAINT_SYSTEM<T,D> )) &, constraint_system ))
        (( typename ARRAY<T>&, constraint_rhs ))
    )
public:
    typedef void result_type;
    void operator()() const
    {
        const int n_embedding = embedding_subsys.linear_index_of_stencil_index.Size();
        embedding_subsys.Init_Stencil_Index_Of_Linear_Index();
        embedding_subsys.stencils.Exact_Resize(n_embedding, false); // uninit'ed
        embedding_subsys.Zero_Stencils();
        const int n_constraint = constraint_system.cell_linear_index_of_stencil_index.Size();
        constraint_system.Init_Stencil_Index_Of_Cell_Linear_Index();
        constraint_system.stencils.Exact_Resize(n_constraint);
        constraint_rhs.Exact_Resize(n_constraint);
    }
};

} // namespace

#define EXPLICIT_INSTANTIATION( T, D ) \
    EXPLICIT_INSTANTIATION_HELPER( T, D, ( DOMAIN_EMBEDDING_CUBE_SUBSYS< T ) ( D > ) )
#define EXPLICIT_INSTANTIATION_HELPER( T, D, T_EMBEDDING_SUBSYS ) \
template int \
Build_Dirichlet_System< T, D, BOOST_PP_SEQ_ENUM( T_EMBEDDING_SUBSYS ) >( \
    EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS const & problem, \
    const MAIN_PARAMS<T,D>& main_params, \
    const ARRAY_VIEW<const T> phi_of_fine_index, \
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys, \
    BOOST_PP_SEQ_ENUM( T_EMBEDDING_SUBSYS )& embedding_subsys, \
    ARRAY_VIEW<T> system_rhs, \
    DIRICHLET_CONSTRAINT_SYSTEM<T,D>& constraint_system, \
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
