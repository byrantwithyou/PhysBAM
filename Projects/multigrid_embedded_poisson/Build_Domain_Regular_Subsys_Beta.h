//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_DOMAIN_REGULAR_SUBSYS_BETA_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_DOMAIN_REGULAR_SUBSYS_BETA_HPP

#include <cassert>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/Functional/APPLY_ASSIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/COMPOSE_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Grid/CELL_VALUE_VIA_AVERAGE_VERTEX_VALUE.h>
#include <Jeffrey_Utilities/Grid/Visit_Cells_With_Sign_Via_Cell_Sign.h>
#include <Jeffrey_Utilities/Grid/VISIT_IF_SIGN_PREDICATE_GRID_VISITOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_X_FUNCTION.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
class DOMAIN_REGULAR_CROSS_SUBSYS;

template< class T, int D >
inline void
Build_Domain_Regular_Subsys_Beta(
    typename EXAMPLE_PARAMS<T,D>::DOMAIN_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys,
    const int domain_sign)
{
    const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
    const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
    const VECTOR<T,D> min_x = As_Vector(main_params.grid.min_x);
    const VECTOR<T,D> max_x = As_Vector(main_params.grid.max_x);
    assert(regular_subsys.beta_of_cell_index.Size() == cell_multi_index_bound.Size());
    Visit_Cells_With_Sign_Via_Cell_Sign_MT(
        main_params.general.n_thread,
        cell_multi_index_bound.Size(),
        As_Const_Array_View(regular_subsys.sign_of_cell_index),
        Make_Visit_If_Sign_Predicate_Grid_Visitor(
            Make_Equal_Function(domain_sign),
            Make_Apply_Assign_Function(
                Make_Array_Wrapper_Function(regular_subsys.beta_of_cell_index),
                Make_Compose_Function(
                    Make_Cell_Value_Via_Average_Vertex_Value(
                        Make_Compose_Function(
                            problem.beta,
                            Make_Multi_Index_X_Function(min_x, max_x, multi_index_bound)
                        )
                    ),
                    cell_multi_index_bound
                )
            )
        )
    );
}

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_DOMAIN_REGULAR_SUBSYS_BETA_HPP
