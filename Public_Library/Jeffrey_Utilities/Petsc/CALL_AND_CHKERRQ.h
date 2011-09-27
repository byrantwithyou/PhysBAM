//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_CALL_AND_CHKERRQ_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_CALL_AND_CHKERRQ_HPP

#ifdef PHYSBAM_USE_PETSC

#define PHYSBAM_PETSC_CALL_AND_CHKERRQ( CallExpr ) \
    do { \
        PetscErrorCode petsc_error_code = CallExpr ; \
        CHKERRQ( petsc_error_code ); \
    } while(false)

#endif // #ifdef PHYSBAM_USE_PETSC

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_CALL_AND_CHKERRQ_HPP
