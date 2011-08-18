//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_DIRICHLET_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_DIRICHLET_SYSTEM_HPP

#include <iosfwd>

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D > class DOMAIN_REGULAR_CROSS_SUBSYS; 
template< class T, int D > struct DIRICHLET_CONSTRAINT_SYSTEM;

template< class T, int D, class T_SIGN, class T_EMBEDDED_SUBSYS >
int Build_Dirichlet_System(
    typename EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    const ARRAY_VIEW<const T> phi_of_fine_index,
    const ARRAY_VIEW<const T_SIGN> sign_of_cell_index,
    DOMAIN_REGULAR_CROSS_SUBSYS<T,D>& regular_subsys,
    T_EMBEDDED_SUBSYS& embedded_subsys,
    ARRAY_VIEW<T> system_rhs,
    DIRICHLET_CONSTRAINT_SYSTEM<T,D>& constraint_system,
    ARRAY<T>& constraint_rhs,
    std::ostream& lout);

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_BUILD_DIRICHLET_SYSTEM_HPP
