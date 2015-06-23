//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_COLLISION_IMPLICIT_OBJECT__
#define __MPM_COLLISION_IMPLICIT_OBJECT__
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{

template<class TV>
class MPM_COLLISION_IMPLICIT_OBJECT:public MPM_COLLISION_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    FRAME<TV>* iot_frame;
    IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* iot;
public:
    using typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE;
    IMPLICIT_OBJECT<TV>* io;
    std::function<FRAME<TV>(T)> func_frame;
    std::function<TWIST<TV>(T)> func_twist;

    MPM_COLLISION_IMPLICIT_OBJECT(IMPLICIT_OBJECT<TV>* io,COLLISION_TYPE type_input,T friction_input);
    virtual ~MPM_COLLISION_IMPLICIT_OBJECT();
    T Phi(const TV& X,T time) const override;
    TV Normal(const TV& X,T time) const override;
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X,T time) const override;
    TV Velocity(const TV& X,T time) const override;
    IMPLICIT_OBJECT<TV>* Get_Implicit_Object(T time) override;
//#####################################################################
};
}
#endif
