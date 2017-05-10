//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __WENO_INTERPOLATION__
#define __WENO_INTERPOLATION__

namespace PhysBAM{

// Interpolate the data at location x.
// If x=0, then interpolates at z2; x=1 interpolates at z3
// Assumes x in [0,1].
// eps is used to prevent division by zero; the paper uses eps=1e-6.
template<class T> T
WENO_Interpolation(T x,T z0,T z1,T z2,T z3,T z4,T z5,T eps);
}
#endif
