//#####################################################################
// Copyright 2015, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_COLLISION_OBJECT__
#define __MPM_COLLISION_OBJECT__
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;
template<class TV> class MPM_COLLISION_OBJECT_PATH;

template<class TV>
class MPM_COLLISION_OBJECT
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:
    IMPLICIT_OBJECT<TV>* impob;
    MPM_COLLISION_OBJECT_PATH<TV>* path;
    const bool sticky;
    const T friction;

    MPM_COLLISION_OBJECT(IMPLICIT_OBJECT<TV>* impob,MPM_COLLISION_OBJECT_PATH<TV>* path,const bool sticky=false,const T friction=0);
    virtual ~MPM_COLLISION_OBJECT();

    bool Detect(const T t,const TV& x,T* phi=0,TV* n=0) const;
    bool Collide(const T t,const TV& x,TV& v,T* phi=0,TV* n=0,bool apply_friction=true) const;
    void Collide_Static(const T t,const TV& x,const TV& n,TV& v,bool apply_friction=true) const;
//#####################################################################
};
}
#endif
