//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_COLLISION_OBJECT__
#define __MPM_COLLISION_OBJECT__
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
namespace PhysBAM{
template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class MPM_COLLISION_OBJECT
{
    typedef typename TV::SCALAR T;
public:
    enum COLLISION_TYPE {stick,slip,separate} type;
    T friction;

    MPM_COLLISION_OBJECT(COLLISION_TYPE type,T friction): type(type),friction(friction) {}
    virtual ~MPM_COLLISION_OBJECT();
    virtual T Phi(const TV& X,T time) const=0;
    virtual TV Normal(const TV& X,T time) const=0;
    virtual SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X,T time) const=0;
    virtual TV Velocity(const TV& X,T time) const=0;
    virtual IMPLICIT_OBJECT<TV>* Get_Implicit_Object(T time)=0;
//#####################################################################
};
}
#endif
