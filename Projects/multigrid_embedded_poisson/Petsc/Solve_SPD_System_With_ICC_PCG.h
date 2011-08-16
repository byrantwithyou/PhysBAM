//#####################################################################
// Copyright 2010, Jeffrey Hellrung, Calvin Wang
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
#define PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP

#include <petsc.h>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>

#include "SYSTEM_REFERENCE.h"

namespace PhysBAM
{

namespace Multigrid_Embedded_Poisson
{

namespace Petsc
{

template< class T >
PetscErrorCode
Solve_SPD_System_With_ICC_PCG(
    const unsigned int n_thread,
    const SYSTEM_REFERENCE<T> system,
    const ARRAY_VIEW<const T> rhs,
    const bool has_constant_vectors_in_null_space,
    unsigned int max_iterations,
    const float relative_tolerance,
    const float absolute_tolerance,
    const bool print_residuals,
    const bool precondition,
    ARRAY_VIEW<T> u_approx);

} // namespace Petsc

} // namespace Multigrid_Embedded_Poisson

} // namespace PhysBAM

#endif // #endif // PHYSBAM_PROJECTS_MULTIGRID_EMBEDDED_POISSON_3D_V2_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
