//#####################################################################
// Copyright 2010, Jeffrey Hellrung, Calvin Wang
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP

#ifdef PHYSBAM_USE_PETSC

#include <iosfwd>

#include <petsc.h>

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <Jeffrey_Utilities/GENERIC_SYSTEM_REFERENCE.h>
#include <Jeffrey_Utilities/ONSTREAM.h>
#include <Jeffrey_Utilities/SOLVER_PARAMS.h>

namespace PhysBAM
{

namespace Petsc
{

template< class T >
PetscErrorCode
Solve_SPD_System_With_ICC_PCG(
    const unsigned int n_thread,
    const SOLVER_PARAMS params,
    const bool has_constant_vectors_in_null_space,
    const GENERIC_SYSTEM_REFERENCE<T> system,
    const ARRAY_VIEW<const T> rhs,
    ARRAY_VIEW<T> x,
    std::ostream& lout = PhysBAM::nout);

} // namespace Petsc

} // namespace PhysBAM

#endif // #ifdef PHYSBAM_USE_PETSC

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SOLVE_SPD_SYSTEM_WITH_ICC_PCG_HPP
