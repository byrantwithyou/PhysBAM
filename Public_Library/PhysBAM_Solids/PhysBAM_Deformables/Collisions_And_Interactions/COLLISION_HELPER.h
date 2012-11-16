//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COLLISION_HELPER__
#define __COLLISION_HELPER__

#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV,class T,int d> TV Compute_Collision_Impulse(const TV& normal,const SYMMETRIC_MATRIX<T,d>& impulse_factor,
    const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction,bool* applied_sticking_impulse);
}
#endif
