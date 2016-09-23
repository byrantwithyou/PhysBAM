//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PHYSBAM_ATTRIBUTE__
#define __PHYSBAM_ATTRIBUTE__

#ifdef _WIN32
#  define PHYSBAM_ALWAYS_INLINE
#  define PHYSBAM_FLATTEN
#else
#  ifdef NDEBUG
#    define PHYSBAM_ALWAYS_INLINE __attribute__ ((always_inline))
#    ifdef __clang__
#      define PHYSBAM_FLATTEN
#    else
#      define PHYSBAM_FLATTEN __attribute__ ((flatten))
#    endif
#  else
#    define PHYSBAM_ALWAYS_INLINE
#    define PHYSBAM_FLATTEN
#  endif
#endif

#endif
