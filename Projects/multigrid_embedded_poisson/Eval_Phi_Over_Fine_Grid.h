//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_EVAL_PHI_OVER_FINE_GRID_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_EVAL_PHI_OVER_FINE_GRID_HPP

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

#include "Params/MAIN_PARAMS.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

template< class T, int D >
void
Eval_Phi_Over_Fine_Grid(
    const MAIN_PARAMS<T,D>& main_params,
    ARRAY_VIEW<T> phi_of_fine_index,
    const int sign_of_zero = -1);

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_V2_EVAL_PHI_OVER_FINE_GRID_HPP
