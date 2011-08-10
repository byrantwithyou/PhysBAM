//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_DOMAIN_REGULAR_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_DOMAIN_REGULAR_SUBSYS_HPP

#include <cassert>

#include <iosfwd>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/Functional/APPLY_ASSIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Grid/CELL_VALUE_VIA_AVERAGE_VERTEX_VALUE.h>
#include <Jeffrey_Utilities/Grid/Visit_Cells_With_Sign_Via_Cell_Sign.h>
#include <Jeffrey_Utilities/Grid/VISIT_IF_SIGN_PREDICATE_GRID_VISITOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/ONSTREAM.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Build_Domain_Regular_Subsys_Rhs.h"
#include "DOMAIN_REGULAR_CROSS_SUBSYS.h"
#include "INIT_CROSS_STENCIL_CELL_VISITOR.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D, class T_BETA_OF_INDEX, class T_F_OF_INDEX >
int
Build_Domain_Regular_Subsys(
    const unsigned int n_thread,
    const VECTOR<T,D> dx,
    const MULTI_INDEX_BOUND<D> cell_multi_index_bound,
    const T_BETA_OF_INDEX& beta_of_index,
    const T_F_OF_INDEX& f_of_index,
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys,
    ARRAY_VIEW<T> system_rhs,
    std::ostream& lout = PhysBAM::nout)
{
    assert(n_thread >= 1);

    BASIC_TIMER timer;

    assert(regular_subsys.sign_of_cell_index.Size() == cell_multi_index_bound.Size());
    assert(regular_subsys.beta_of_cell_index.Size() == cell_multi_index_bound.Size());
    assert(regular_subsys.stencil_of_index.Size() == (cell_multi_index_bound + 1).Size());
    assert(system_rhs.Size() == (cell_multi_index_bound + 1).Size());

    lout << "Initializing beta and regular subsystem stencils...";
    lout.flush();
    timer.Restart();
    regular_subsys.Zero_Stencils(n_thread);
    Visit_Cells_With_Sign_Via_Cell_Sign_MT(
        n_thread,
        cell_multi_index_bound.Size(),
        As_Const_Array_View(regular_subsys.sign_of_cell_index),
        Make_Visit_If_Sign_Predicate_Grid_Visitor(
            Make_Equal_Function(-1),
            Make_Apply_Assign_Function(
                Make_Array_Wrapper_Function(regular_subsys.beta_of_cell_index),
                Make_Compose_Function(
                    Make_Cell_Value_Via_Average_Vertex_Value(beta_of_index),
                    cell_multi_index_bound
                )
            )
        )
    );
    // TODO: MT
    Visit_Cells_With_Sign_Via_Cell_Sign(
        cell_multi_index_bound.Size(),
        As_Const_Array_View(regular_subsys.sign_of_cell_index),
        Make_Visit_If_Sign_Predicate_Grid_Visitor(
            Make_Equal_Function(-1),
            Make_Init_Cross_Stencil_Cell_Visitor(
                dx, cell_multi_index_bound,
                As_Const_Array_View(regular_subsys.beta_of_cell_index),
                Make_Array_Wrapper_Function(regular_subsys.stencil_of_index)
            )
        )
    );
    lout << timer.Elapsed() << " s" << std::endl;

    lout << "Initializing regular subsystem rhs...";
    lout.flush();
    timer.Restart();
    Build_Domain_Regular_Subsys_Rhs(
        n_thread,
        dx.Product(), cell_multi_index_bound,
        -1, // domain_sign
        As_Const_Array_View(regular_subsys.sign_of_cell_index),
        Make_Compose_Function(
            Make_Cell_Value_Via_Average_Vertex_Value(f_of_index),
            cell_multi_index_bound
        ),
        system_rhs
    );
    lout << timer.Elapsed() << " s" << std::endl;

    return 0;
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_DOMAIN_REGULAR_SUBSYS_HPP
