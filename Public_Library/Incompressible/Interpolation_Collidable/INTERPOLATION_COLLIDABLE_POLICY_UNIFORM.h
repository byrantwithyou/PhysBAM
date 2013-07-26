//#####################################################################
// Copyright 2009, Andrew Selle, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERPOLATION_COLLIDABLE_POLICY_UNIFORM__
#define __INTERPOLATION_COLLIDABLE_POLICY_UNIFORM__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <Incompressible/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <Incompressible/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD.h>
namespace PhysBAM{

template<class TV>
struct INTERPOLATION_COLLIDABLE_POLICY
{
private:
    typedef typename TV::SCALAR T;
public:
    typedef AVERAGING_UNIFORM<TV> AVERAGING;
    typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> FACE_LOOKUP_COLLIDABLE;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<TV,T> LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR;
    typedef LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM<TV,T> LINEAR_INTERPOLATION_COLLIDABLE_FACE_SCALAR;
    // slip collidable
    typedef FACE_LOOKUP_UNIFORM<TV> FACE_LOOKUP_SLIP;
    typedef FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<TV> FACE_LOOKUP_COLLIDABLE_SLIP;
};
//#####################################################################
}
#endif
