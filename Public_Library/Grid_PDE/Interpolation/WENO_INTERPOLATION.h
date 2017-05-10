//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __WENO_INTERPOLATION__
#define __WENO_INTERPOLATION__

#include <Core/Arrays_Nd/ARRAYS_ND_VIEW.h>
namespace PhysBAM{

template<class T,class TV_INT> class ARRAYS_ND_VIEW;

// Interpolate the data at location x.
// If x=0, then interpolates at z2; x=1 interpolates at z3
// Assumes x in [0,1].
// eps is used to prevent division by zero; the paper uses eps=1e-6.
template<class T,class U> U
WENO_Interpolation(T x,U z0,U z1,U z2,U z3,U z4,U z5,T eps);

template<class T,class U> inline U
WENO_Interpolation(T x,U* z,T eps)
{
    return WENO_Interpolation(x,z[0],z[1],z[2],z[3],z[4],z[5],eps);
}

template<class T,class U,int d> U
WENO_Interpolation(const VECTOR<T,d>& x,const ARRAYS_ND_BASE<U,VECTOR<int,d> >& z,const VECTOR<int,d>& index,T eps);
}
#endif
