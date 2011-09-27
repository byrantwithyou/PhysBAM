//#####################################################################
// Copyright 2010, Jeffrey Hellrung, Calvin Wang
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_PRINT_KSP_INFO_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_PRINT_KSP_INFO_HPP

#ifdef PHYSBAM_USE_PETSC

#include <iosfwd>

#include <petsc.h>
#include <petscksp.h>

#include <Jeffrey_Utilities/ONSTREAM.h>

namespace PhysBAM
{

namespace Petsc
{

PetscErrorCode
Print_KSP_Info(
    const KSP petsc_ksp,
    std::ostream& lout = PhysBAM::nout);

} // namespace Petsc

} // namespace PhysBAM

#endif // #ifdef PHYSBAM_USE_PETSC

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_PRINT_KSP_INFO_HPP
