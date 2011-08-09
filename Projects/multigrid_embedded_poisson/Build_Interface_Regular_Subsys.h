//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_INTERFACE_REGULAR_SUBSYS_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_INTERFACE_REGULAR_SUBSYS_HPP

#include <cassert>

#include <iostream>

#include <Jeffrey_Utilities/ARRAY_OPS.h>
#include <Jeffrey_Utilities/BASIC_TIMER.h>
#include <Jeffrey_Utilities/Functional/APPLY_ASSIGN_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/ARRAY_WRAPPER_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Functional/NOT_EQUAL_FUNCTION.h>
#include <Jeffrey_Utilities/Grid/Visit_Cells_With_Sign_Via_Cell_Sign.h>
#include <Jeffrey_Utilities/Grid/VISIT_IF_SIGN_PREDICATE_GRID_VISITOR.h>
#include <Jeffrey_Utilities/Multi_Index/MULTI_INDEX_BOUND.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

#include "Build_Domain_Regular_Subsys_Beta.h"
#include "Build_Domain_Regular_Subsys_Rhs.h"
#include "INIT_CROSS_CONSTBETA_STENCIL_CELL_VISITOR.h"
#include "INIT_CROSS_STENCIL_CELL_VISITOR.h"
#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D > class DOMAIN_REGULAR_CROSS_SUBSYS;

namespace Detail_Build_Interface_Regular_Subsys
{

template< class T, int D >
inline void
Init_Beta_And_Stencils(
    typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys);

} // namespace Detail_Build_Interface_Regular_Subsys

template< class T, int D, class T_REGULAR_SUBSYS >
int
Build_Interface_Regular_Subsys(
    typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    const ARRAY_VIEW<const T> phi_of_fine_index,
    T_REGULAR_SUBSYS& regular_subsys,
    ARRAY_VIEW<T> system_rhs)
{
    BASIC_TIMER timer;

    const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
    const MULTI_INDEX_BOUND<D> multi_index_bound = cell_multi_index_bound + 1;
    const MULTI_INDEX_BOUND<D> fine_multi_index_bound = 2 * multi_index_bound - 1;

    assert(phi_of_fine_index.Size() == fine_multi_index_bound.Size());
    assert(regular_subsys.sign_of_cell_index.Size() == cell_multi_index_bound.Size());
    assert(regular_subsys.stencil_of_index.Size() >= multi_index_bound.Size());
    assert(system_rhs.Size() >= multi_index_bound.Size());

    std::cout << "Initializing beta and regular subsystem stencils...";
    std::cout.flush();
    timer.Restart();
    regular_subsys.Zero_Stencils_MT(main_params.general.n_thread);
    Detail_Build_Interface_Regular_Subsys::Init_Beta_And_Stencils(problem, main_params, regular_subsys);
    std::cout << timer.Elapsed() << " s" << std::endl;

    std::cout << "Initializing regular subsystem rhs...";
    std::cout.flush();
    timer.Restart();
    Build_Domain_Regular_Subsys_Rhs(problem.negative, main_params, regular_subsys, system_rhs, -1);
    Build_Domain_Regular_Subsys_Rhs(problem.positive, main_params, regular_subsys, system_rhs, +1);
    std::cout << timer.Elapsed() << " s" << std::endl;

    return 0;
}

namespace Detail_Build_Interface_Regular_Subsys
{

template< class T, int D >
inline void
Init_Beta_And_Stencils(
    typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys)
{
    const MULTI_INDEX_BOUND<D> cell_multi_index_bound(As_Vector<int>(main_params.grid.n_cell));
    assert(regular_subsys.beta_of_cell_index.Size() == cell_multi_index_bound.Size());

    Build_Domain_Regular_Subsys_Beta(problem.negative, main_params, regular_subsys, -1);
    Build_Domain_Regular_Subsys_Beta(problem.positive, main_params, regular_subsys, +1);

    // TODO: MT
    Visit_Cells_With_Sign_Via_Cell_Sign(
        cell_multi_index_bound.Size(),
        As_Const_Array_View(regular_subsys.sign_of_cell_index),
        Make_Visit_If_Sign_Predicate_Grid_Visitor(
            Make_Not_Equal_Function(0),
            Make_Init_Cross_Stencil_Cell_Visitor(
                regular_subsys.dx,
                cell_multi_index_bound,
                As_Const_Array_View(regular_subsys.beta_of_cell_index),
                Make_Array_Wrapper_Function(regular_subsys.stencil_of_index)
            )
        )
    );
}

} // namespace Detail_Build_Interface_Regular_Subsys

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_INTERFACE_REGULAR_SUBSYS_HPP
