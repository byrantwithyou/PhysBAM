//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __constants__
#define __constants__

#include <cmath>
namespace PhysBAM{

static const constexpr double pi=M_PI;

static const constexpr double speed_of_light=2.99792458e8; // m/s
static const constexpr double plancks_constant=6.6260755e-34; // J*s
static const constexpr double boltzmanns_constant=1.380658e-23; // J/K

static const constexpr double unit_sphere_size[4]={0,2,pi,pi*4/3};

}
#endif
