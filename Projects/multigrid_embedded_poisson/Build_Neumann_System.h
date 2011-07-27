//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_NEUMANN_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_NEUMANN_SYSTEM_HPP

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D, class T_REGULAR_SUBSYS, class T_EMBEDDED_SUBSYS >
int Build_Neumann_System(
    typename EXAMPLE_PARAMS<T,D>::NEUMANN_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    const ARRAY_VIEW<const T> phi_of_fine_index,
    T_REGULAR_SUBSYS& regular_subsys,
    T_EMBEDDED_SUBSYS& embedded_subsys,
    ARRAY_VIEW<T> system_rhs);

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_NEUMANN_SYSTEM_HPP
