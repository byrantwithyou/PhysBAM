//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_COUPLED_IMPULSE
//#####################################################################
#ifndef __COMBINED_COLLISIONS_COUPLED_IMPULSE__
#define __COMBINED_COLLISIONS_COUPLED_IMPULSE__
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_DEFORMABLE_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/COMBINED_COLLISIONS_RIGID_IMPULSE.h>

namespace PhysBAM{

template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class SOLID_BODY_COLLECTION;

template<class TV>
struct COMBINED_COLLISIONS_COUPLED_IMPULSE:public COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>,public COMBINED_COLLISIONS_RIGID_IMPULSE<TV>
{
    typedef typename TV::SCALAR T;
    typedef COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV> DBASE;using DBASE::impulse;
    typedef COMBINED_COLLISIONS_RIGID_IMPULSE<TV> RBASE;using RBASE::wrench;
    using RBASE::Apply_Rigid_Impulse;

    COMBINED_COLLISIONS_COUPLED_IMPULSE(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,bool do_collision);
    virtual ~COMBINED_COLLISIONS_COUPLED_IMPULSE();

    virtual void Resize();
    virtual void Apply(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time);
    virtual void Clear(const ARRAY<COMBINED_BODY_ID>& list);
    virtual void Clear();
    virtual void Setup_Discover_State(T dt,T time);
    virtual void Setup_Impulse_State(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time);
    virtual void Finish_Step(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time);
    virtual T Kinetic_Energy(const ARRAY<COMBINED_BODY_ID>& list) const;

    void Apply_Coupled_Impulse(int p,int b,const TV& j);
};
}
#endif
