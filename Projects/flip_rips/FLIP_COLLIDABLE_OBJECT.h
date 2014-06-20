//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLIP_COLLIDABLE_OBJECT__
#define __FLIP_COLLIDABLE_OBJECT__

#include <Tools/Vectors/VECTOR.h>

#include "FLIP_TRANSLATION_PATH.h"
#include "FLIP_ROTATION_PATH.h"
#include "FLIP_LEVELSET.h"

namespace PhysBAM{

template<class TV>
class FLIP_COLLIDABLE_OBJECT
{
    typedef typename TV::SCALAR T;

public:

    IMPLICIT_OBJECT<TV>* ls;
    FLIP_TRANSLATION_PATH<TV>* tp;
    FLIP_ROTATION_PATH<TV>* rp;

    FLIP_COLLIDABLE_OBJECT(
        IMPLICIT_OBJECT<TV>* ls_input,
        FLIP_TRANSLATION_PATH<TV>* tp_input,
        FLIP_ROTATION_PATH<TV>* rp_input)
        :ls(ls_input),tp(tp_input),
         rp(rp_input)
    {}

    virtual ~FLIP_COLLIDABLE_OBJECT()
    {
        delete ls;
        delete tp;
        delete rp;
    }

    bool Collide(const TV& x,const T t,TV& v) const
    {
        TV normal;
        if(Detect_Collision(x,t,&normal)){
            TV V=Velocity(x,t);
            v-=V;
            Collide_Static(x,t,normal,v);
            v+=V;
            return true;}
        return false;
    }

    bool Detect_Collision(const TV& x,const T t,TV* normal=0) const
    {
        TV X=rp->R(t).Inverse_Rotate(x-tp->X(t));
        if(ls->Extended_Phi(X)<0){
            if(normal) *normal=rp->R(t).Rotate(ls->Normal(X));
            return true;}
        return false;
    }

    virtual TV Velocity(const TV& x,const T t) const
    {return TV::Cross_Product(rp->W(t),x-tp->X(t))+tp->V(t);}

protected:

    void Collide_Static(const TV& x,const T t,const TV& normal,TV& v) const
    {
        const T projection=TV::Dot_Product(normal,v);
        if(projection<0) v-=normal*projection;
    }
};

template<class TV>
class FLIP_COLLIDABLE_OBJECT_STATIC: public FLIP_COLLIDABLE_OBJECT<TV>
{
    typedef FLIP_COLLIDABLE_OBJECT<TV> BASE;
    typedef typename TV::SCALAR T;
    
public:

    FLIP_COLLIDABLE_OBJECT_STATIC(IMPLICIT_OBJECT<TV>* const ls_input,const TV& tp_input,const ROTATION<TV>& rp_input)
        :BASE(ls_input,new FLIP_TRANSLATION_PATH_STATIC<TV>(tp_input),new FLIP_ROTATION_PATH_STATIC<TV>(rp_input))
    {}

    virtual ~FLIP_COLLIDABLE_OBJECT_STATIC(){}  
};
}
#endif
