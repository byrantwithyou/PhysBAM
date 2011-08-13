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
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <Jeffrey_Utilities/VISITOR_SEQUENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "BETA_GRAD_U_DOT_N.h"
#include "Build_Interface_Embedding_Subsys.h"
#include "Build_Interface_Regular_Subsys.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "EMBEDDING_UNSTRUCTURED_SUBSYS.h"
#include "Init_Cell_Local_Embedding_Interface_System.h"
#include "INTERFACE_CONSTRAINT_SYSTEM.h"
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

template< class T, int D >
int Build_Interface_System(
    const typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS& problem,
    const MAIN_PARAMS<T,D>& main_params,
    const ARRAY_VIEW<const T> phi_of_fine_index,
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys,
    EMBEDDING_UNSTRUCTURED_SUBSYS<T>& embedding_subsys,
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
        EMBEDDING_UNSTRUCTURED_SUBSYS<T>&
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

    Build_Interface_Embedding_Subsys(
        main_params.general.n_thread,
        min_x, max_x, multi_index_bound,
        As_Const_Array_View(regular_subsys.sign_of_cell_index),
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
        embedding_subsys, system_rhs,
        constraint_system, constraint_rhs,
        lout
    );

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

#define EXPLICIT_INSTANTIATION( T, D ) \
template int \
Build_Interface_System< T, D >( \
    EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS const & problem, \
    const MAIN_PARAMS<T,D>& main_params, \
    const ARRAY_VIEW<const T> phi_of_fine_index, \
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys, \
    EMBEDDING_UNSTRUCTURED_SUBSYS<T>& embedding_subsys, \
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
