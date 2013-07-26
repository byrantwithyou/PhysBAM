//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FIRE_INTERPOLATION_POLICY__
#define __FIRE_INTERPOLATION_POLICY__

#include <Dynamics/Interpolation/FIRE_INTERPOLATION_FORWARD.h>
namespace PhysBAM{

template<class T,int d> class VECTOR;

template<class TV>
struct FIRE_INTERPOLATION_POLICY
{
private:
    typedef typename TV::SCALAR T;
public:
    // multiphase fire
    typedef FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<TV> FACE_LOOKUP_FIRE_MULTIPHASE;
    typedef AVERAGING_UNIFORM<TV,FACE_LOOKUP_FIRE_MULTIPHASE> AVERAGING_FIRE_MULTIPHASE;
    typedef LINEAR_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_FIRE_MULTIPHASE> LINEAR_INTERPOLATION_FIRE_MULTIPHASE;
    // multiphase fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<TV,FACE_LOOKUP_FIRE_MULTIPHASE> FACE_LOOKUP_FIRE_MULTIPHASE_COLLIDABLE;
};
//#####################################################################
}
#endif
