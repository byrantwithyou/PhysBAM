//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_RIGID_IMPULSE
//#####################################################################
#ifndef __COMBINED_COLLISIONS_RIGID_IMPULSE__
#define __COMBINED_COLLISIONS_RIGID_IMPULSE__
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>

namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLISIONS;
template<class TV> class RIGID_BODY_COLLECTION;

template<class TV>
struct COMBINED_COLLISIONS_RIGID_IMPULSE:public virtual COMBINED_COLLISIONS<TV>::IMPULSE
{
    typedef typename TV::SCALAR T;
    ARRAY<TWIST<TV> > wrench;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions;
    bool update_positions_for_collisions_on_apply;
    bool update_positions_for_contact_on_apply;
    bool use_collisions;

    COMBINED_COLLISIONS_RIGID_IMPULSE(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions_input,bool do_collision);
    virtual ~COMBINED_COLLISIONS_RIGID_IMPULSE();

    virtual void Resize();
    virtual void Apply(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time);
    virtual void Clear(const ARRAY<COMBINED_BODY_ID>& list);
    virtual void Clear();
    virtual void Setup_Discover_State(T dt,T time);
    virtual void Setup_Impulse_State(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time);
    virtual void Finish_Step(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time);
    virtual T Kinetic_Energy(const ARRAY<COMBINED_BODY_ID>& list) const;

    void Apply_Rigid_Impulse(int b,const TV& location,const TWIST<TV>& j);
    T Kinetic_Energy_Change(int b,const TV& location,const TWIST<TV>& j) const;
    void Apply_Rigid_Impulse(int a,int b,const TV& location,const TWIST<TV>& j);
    TWIST<TV> Body_Twist(int b) const;
    TV Velocity_At_Point(int e,const TV& location) const;
    void Apply(int b,T dt,T time);
    void Setup_Impulse_State(int b,T dt,T time);
    void Finish_Step(int b,T dt,T time);
    T Kinetic_Energy(int b) const;
};
}
#endif
