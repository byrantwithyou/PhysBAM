//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SCOPED_FINALIZE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SCOPED_FINALIZE_HPP

#ifdef PHYSBAM_USE_PETSC

#include <boost/preprocessor/cat.hpp>

#ifdef __COUNTER__

#define PHYSBAM_PETSC_SCOPED_FINALIZE() \
    PHYSBAM_PETSC_SCOPED_FINALIZE_impl( BOOST_PP_CAT( _petsc_finalize_, __COUNTER__ ) )
#define PHYSBAM_PETSC_SCOPED_FINALIZE_IF( Cond ) \
    PHYSBAM_PETSC_SCOPED_FINALIZE_IF_impl( Cond, BOOST_PP_CAT( _petsc_finalize_, __COUNTER__ ) )

#else // #ifdef __COUNTER__

#define PHYSBAM_PETSC_SCOPED_FINALIZE() \
    PHYSBAM_PETSC_SCOPED_FINALIZE_impl( BOOST_PP_CAT( _petsc_finalize_, __LINE__ ) )
#define PHYSBAM_PETSC_SCOPED_FINALIZE_IF( Cond ) \
    PHYSBAM_PETSC_SCOPED_FINALIZE_IF_impl( Cond, BOOST_PP_CAT( _petsc_finalize_, __LINE__ ) )

#endif // #ifdef __COUNTER__

#define PHYSBAM_PETSC_SCOPED_FINALIZE_impl( RAIIName ) \
    struct BOOST_PP_CAT( RAIIName, _t ) { \
        ~ BOOST_PP_CAT( RAIIName, _t ) () { PetscFinalize(); } \
    } RAIIName; static_cast<void>( RAIIName )

#define PHYSBAM_PETSC_SCOPED_FINALIZE_IF_impl( Cond, RAIIName ) \
    struct BOOST_PP_CAT( RAIIName, _t ) { \
        bool b; \
        ~ BOOST_PP_CAT( RAIIName, _t ) () { if(b) PetscFinalize(); } \
    } RAIIName = { Cond }; static_cast<void>( RAIIName )

#endif // #ifdef PHYSBAM_USE_PETSC

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SCOPED_FINALIZE_HPP
