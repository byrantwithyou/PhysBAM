//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_AND_SOLVE_DIRICHLET_SYSTEM_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_AND_SOLVE_DIRICHLET_SYSTEM_HPP

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

#include "Params/EXAMPLE_PARAMS.h"
#include "Params/MAIN_PARAMS.h"
#include "RAND_MT19937_UNIFORM_REAL.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

template< class T, int D >
int Build_And_Solve_Dirichlet_System(
    typename EXAMPLE_PARAMS<T,D>::DIRICHLET_PARAMS const & problem,
    const MAIN_PARAMS<T,D>& main_params,
    typename RAND_MT19937_UNIFORM_REAL<T>::type& rand,
    const ARRAY_VIEW<const T> phi_of_fine_index);

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_BUILD_AND_SOLVE_DIRICHLET_SYSTEM_HPP
