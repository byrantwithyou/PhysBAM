//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_BUILD_AND_SOLVE_INTERFACE_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_BUILD_AND_SOLVE_INTERFACE_SYSTEM_HPP

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "RAND_MT19937_UNIFORM_REAL.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D >
int Build_And_Solve_Interface_System(
    const typename EXAMPLE_PARAMS<T,D>::INTERFACE_PARAMS& problem,
    const MAIN_PARAMS<T,D>& main_params,
    typename RAND_MT19937_UNIFORM_REAL<T>::type& rand,
    const ARRAY_VIEW<const T> phi_of_fine_index);

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_BUILD_AND_SOLVE_INTERFACE_SYSTEM_HPP
