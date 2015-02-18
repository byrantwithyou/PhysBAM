//#####################################################################
// Copyright 2015, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_COLLISION_OBJECT_UNIFORM_PATH__
#define __MPM_COLLISION_OBJECT_UNIFORM_PATH__
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT_PATH.h>
namespace PhysBAM{

template<class TV> class ROTATION;
template<class TV> class FRAME;
template<class TV> class TWIST;

template<class TV>
class MPM_COLLISION_OBJECT_UNIFORM_PATH:public MPM_COLLISION_OBJECT_PATH<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:
    FRAME<TV> frame;
    TWIST<TV> twist;

    MPM_COLLISION_OBJECT_UNIFORM_PATH(const FRAME<TV>& frame,const TWIST<TV>& twist);
    virtual ~MPM_COLLISION_OBJECT_UNIFORM_PATH();

    TV Translation(const T t) const;
    TV Translation_Velocity(const T t) const;
    ROTATION<TV> Rotation(const T t) const;
    T_SPIN Angular_Velocity(const T t) const;
//#####################################################################
};
}
#endif
