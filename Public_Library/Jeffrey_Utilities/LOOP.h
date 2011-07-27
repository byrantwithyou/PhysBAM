//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LOOP_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LOOP_HPP

#include <boost/preprocessor/cat.hpp>

#ifdef __COUNTER__

#define PHYSBAM_LOOP( N ) PHYSBAM_LOOP_impl( N, BOOST_PP_CAT( _physbam_loop_, __COUNTER__ ) )

#else // #ifdef __COUNTER__

#define PHYSBAM_LOOP( N ) PHYSBAM_LOOP_impl( N, BOOST_PP_CAT( _physbam_loop_, __LINE__ ) )

#endif // #ifdef __COUNTER__

#define PHYSBAM_LOOP_impl( N, x ) for(int x = 1; x <= ( N ); ++x)

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_LOOP_HPP
