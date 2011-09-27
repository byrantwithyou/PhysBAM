//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SCOPED_DESTROY_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SCOPED_DESTROY_HPP

#ifdef PHYSBAM_USE_PETSC

#include <boost/preprocessor/cat.hpp>

#ifdef __COUNTER__

#define PHYSBAM_PETSC_SCOPED_DESTROY( Type, X ) \
    PHYSBAM_PETSC_SCOPED_DESTROY_impl( \
        Type, \
        X, \
        BOOST_PP_CAT( _petsc_, BOOST_PP_CAT( Type, BOOST_PP_CAT( _destroy_, __COUNTER__ ) ) ) \
    )
#define PHYSBAM_PETSC_SCOPED_DESTROY_IF( Type, X, Cond ) \
    PHYSBAM_PETSC_SCOPED_DESTROY_IF_impl( \
        Type, \
        X, \
        Cond, \
        BOOST_PP_CAT( _petsc_, BOOST_PP_CAT( Type, BOOST_PP_CAT( _destroy_, __COUNTER__ ) ) ) \
    )

#else // #ifdef __COUNTER__

#define PHYSBAM_PETSC_SCOPED_DESTROY( Type, X ) \
    PHYSBAM_PETSC_SCOPED_DESTROY_impl( \
        Type, \
        X, \
        BOOST_PP_CAT( _petsc_, BOOST_PP_CAT( Type, BOOST_PP_CAT( _destroy_, __LINE__ ) ) ) \
    )
#define PHYSBAM_PETSC_SCOPED_DESTROY_IF( Type, X, Cond ) \
    PHYSBAM_PETSC_SCOPED_DESTROY_IF_impl( \
        Type, \
        X, \
        Cond, \
        BOOST_PP_CAT( _petsc_, BOOST_PP_CAT( Type, BOOST_PP_CAT( _destroy_, __LINE__ ) ) ) \
    )

#endif // #ifdef PHYSBAM_USE_PETSC

#endif // #ifdef __COUNTER__

#define PHYSBAM_PETSC_SCOPED_DESTROY_impl( Type, X, RAIIName ) \
    struct BOOST_PP_CAT( RAIIName, _t ) { \
        const Type x; \
        ~ BOOST_PP_CAT( RAIIName, _t ) () { BOOST_PP_CAT( Type, Destroy ) (x); } \
    } RAIIName = { X }; static_cast<void>( RAIIName )

#define PHYSBAM_PETSC_SCOPED_DESTROY_IF_impl( Type, X, Cond, RAIIName ) \
    struct BOOST_PP_CAT( RAIIName, _t ) { \
        const Type x; bool b; \
        ~ BOOST_PP_CAT( RAIIName, _t ) () { if(b) BOOST_PP_CAT( Type, Destroy ) (x); } \
    } RAIIName = { X, Cond }; static_cast<void>( RAIIName )

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_PETSC_SCOPED_DESTROY_HPP
