//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_COLLISION_IMPLICIT_SPHERE__
#define __MPM_COLLISION_IMPLICIT_SPHERE__
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{

template<class TV>
class MPM_COLLISION_IMPLICIT_SPHERE:public MPM_COLLISION_IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef MPM_COLLISION_IMPLICIT_OBJECT<TV> BASE;
    mutable ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> > *sphere;
    T radius;
    T shrinking_speed;
public:
    using BASE::io;using BASE::func_frame; using BASE::func_twist;
    using typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE;

    MPM_COLLISION_IMPLICIT_SPHERE(COLLISION_TYPE type_input,T friction_input,
        std::function<FRAME<TV>(T)> func_frame=0,std::function<TWIST<TV>(T)> func_twist=0,TV center=TV(),T radius=1,
        T shrinking_speed=0);
    virtual ~MPM_COLLISION_IMPLICIT_SPHERE();
    T Phi(const TV& X,T time) const override;
    TV Normal(const TV& X,T time) const override;
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X,T time) const override;
    TV Velocity(const TV& X,T time) const override;
    IMPLICIT_OBJECT<TV>* Get_Implicit_Object(T time) override;
    void Update(T time) const;
//#####################################################################
};
}
#endif
