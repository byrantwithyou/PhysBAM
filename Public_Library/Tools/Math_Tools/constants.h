//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __constants__
#define __constants__

#include <cmath>
namespace PhysBAM{

const double pi=4*atan(1.);

const double speed_of_light=2.99792458e8; // m/s
const double plancks_constant=6.6260755e-34; // J*s
const double boltzmanns_constant=1.380658e-23; // J/K

const double unit_sphere_size[4]={0,2,pi,pi*4/3};

}
#endif
