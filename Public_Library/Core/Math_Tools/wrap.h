//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Igor Neverov, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function wrap
//#####################################################################
//
// wrap(i,n) adds a multiple of n to i to bring it into the set [0;n)
// i + k*n is in [0,n)
//               
//#####################################################################
#ifndef __wrap__
#define __wrap__

#include <Core/Vectors/VECTOR_FORWARD.h>
#include <cmath>
namespace PhysBAM{

inline int wrap(const int i,const int n)
{
    int k=i%n;
    if(k<0) k+=n;
    return k;
}

// ensure lower <= value < upper
inline int wrap(const int value,const int lower,const int upper)
{return wrap(value-lower,upper-lower)+lower;}

inline float wrap(const float value,const float lower,const float upper)
{float r=fmod(value-lower,upper-lower);if(r<0) return r+upper;else return r+lower;}

inline double wrap(const double value,const double lower,const double upper)
{double r=fmod(value-lower,upper-lower);if(r<0) return r+upper;else return r+lower;}

}
#endif
