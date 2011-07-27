//#####################################################################
// Copyright 2010, Jeffrey Hellrung, Calvin Wang
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP

#include <petsc.h>

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>

#include "../Params/EXAMPLE_PARAMS.h"
#include "../Params/MAIN_PARAMS.h"

namespace PhysBAM
{

namespace Embedded_Poisson_V2
{

namespace Petsc
{

template< class T, int D, class T_SYSTEM >
PetscErrorCode
Solve_SPD_System_With_ICC_PCG(
    const MAIN_PARAMS<T,D>& main_params,
    const T_SYSTEM& system, const ARRAY_VIEW<const T> rhs,
    const bool has_constant_vectors_in_null_space,
    ARRAY_VIEW<T> u_approx);

} // namespace Petsc

} // namespace Embedded_Poisson_V2

} // namespace PhysBAM

#endif // #endif // PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
