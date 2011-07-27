//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_DEFORMABLE_IMPULSE
//#####################################################################
#ifndef __COMBINED_COLLISIONS_DEFORMABLE_IMPULSE__
#define __COMBINED_COLLISIONS_DEFORMABLE_IMPULSE__
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>

namespace PhysBAM{

template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class SOLID_BODY_COLLECTION;

template<class TV>
struct COMBINED_COLLISIONS_DEFORMABLE_IMPULSE:public virtual COMBINED_COLLISIONS<TV>::IMPULSE
{
    typedef typename TV::SCALAR T;typedef typename COMBINED_COLLISIONS<TV>::IMPULSE BASE;

    ARRAY<TV> impulse;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;

    COMBINED_COLLISIONS_DEFORMABLE_IMPULSE(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input);
    virtual ~COMBINED_COLLISIONS_DEFORMABLE_IMPULSE();

    virtual void Resize();
    virtual void Apply(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time);
    virtual void Clear(const ARRAY<COMBINED_BODY_ID>& list);
    virtual void Clear();
    virtual void Setup_Discover_State(T dt,T time);
    virtual void Setup_Impulse_State(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time);
    virtual void Finish_Step(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time);
    virtual T Kinetic_Energy(const ARRAY<COMBINED_BODY_ID>& list) const;

    void Apply(int e,T dt,T time);
    TV Particle_Velocity(int e) const;
    void Setup_Impulse_State(int p,T dt,T time){}
    void Finish_Step(int p,T dt,T time){}
    T Kinetic_Energy(int p) const;
    T Kinetic_Energy_Change(int b,const TV& j) const;
};
}
#endif
