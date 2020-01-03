//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_COLLISION_OBJECT__
#define __MPM_COLLISION_OBJECT__
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{
template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class MPM_COLLISION_OBJECT
{
    typedef typename TV::SCALAR T;
    FRAME<TV>* iot_frame=0;
    IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* iot=0;
public:
    enum COLLISION_TYPE {stick,slip,separate} type;
    T friction;
    IMPLICIT_OBJECT<TV>* io;

    std::function<FRAME<TV>(T)> func_frame;
    std::function<TWIST<TV>(T)> func_twist;

    MPM_COLLISION_OBJECT(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type,T friction):
        type(type),friction(friction),io(io)
    {}

    ~MPM_COLLISION_OBJECT();

    T Phi(const TV& X,T time) const;
    TV Normal(const TV& X,T time) const;
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X,T time) const;
    TV Velocity(const TV& X,T time) const;
    IMPLICIT_OBJECT<TV>* Get_Implicit_Object(T time);
//#####################################################################
};
}
#endif
