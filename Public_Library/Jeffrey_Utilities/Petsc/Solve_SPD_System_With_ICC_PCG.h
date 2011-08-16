//#####################################################################
// Copyright 2010, Jeffrey Hellrung, Calvin Wang
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP

#include <iosfwd>

#include <petsc.h>

#include <Jeffrey_Utilities/ONSTREAM.h>
#include <Jeffrey_Utilities/Petsc/GENERIC_SYSTEM_REFERENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>

namespace PhysBAM
{

namespace Petsc
{

template< class T >
PetscErrorCode
Solve_SPD_System_With_ICC_PCG(
    const unsigned int n_thread,
    const GENERIC_SYSTEM_REFERENCE<T> system,
    const ARRAY_VIEW<const T> rhs,
    const bool has_constant_vectors_in_null_space,
    unsigned int max_iterations,
    const float relative_tolerance,
    const float absolute_tolerance,
    const bool print_residuals,
    const bool precondition,
    ARRAY_VIEW<T> u_approx,
    std::ostream& lout = PhysBAM::nout);

} // namespace Petsc

} // namespace PhysBAM

#endif // #endif // PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
