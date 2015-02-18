//#####################################################################
// Copyright 2015, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_COLLISION_OBJECT_PATH__
#define __MPM_COLLISION_OBJECT_PATH__
namespace PhysBAM{

template<class TV> class ROTATION;

template<class TV>
class MPM_COLLISION_OBJECT_PATH
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:

    virtual ~MPM_COLLISION_OBJECT_PATH();

    virtual TV Translation(const T t) const=0;
    virtual TV Translation_Velocity(const T t) const=0;
    virtual ROTATION<TV> Rotation(const T t) const=0;
    virtual T_SPIN Angular_Velocity(const T t) const=0;

//#####################################################################
};
}
#endif
