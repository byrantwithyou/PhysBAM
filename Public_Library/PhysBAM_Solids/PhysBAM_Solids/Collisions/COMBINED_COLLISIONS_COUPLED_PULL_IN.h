//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_COUPLED_PULL_IN
//#####################################################################
#ifndef __COMBINED_COLLISIONS_COUPLED_PULL_IN__
#define __COMBINED_COLLISIONS_COUPLED_PULL_IN__
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_DEFORMABLE_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/COMBINED_COLLISIONS_RIGID.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/COMBINED_COLLISIONS_COUPLED_IMPULSE.h>

namespace PhysBAM{

template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class COMBINED_COLLISIONS_COUPLED_IMPULSE;
template<class TV> class RIGID_DEFORMABLE_COLLISIONS;

// TODO: Update spatial partition.

template<class TV>
class COMBINED_COLLISIONS_COUPLED_PULL_IN:public virtual COMBINED_COLLISIONS<TV>::COLLIDER
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename COMBINED_COLLISIONS<TV>::IMPULSE IMPULSE;
public:
    RIGID_DEFORMABLE_COLLISIONS<TV>& rigid_deformable_collisions;

    struct COLLISION
    {
        int body;
        int p;
        T target_relative_speed;
        T relative_speed;
        TV normal;
        mutable TV total_impulse;
    };

    ARRAY<COLLISION> collisions;
    T collision_body_thickness;
    const T dt;
    mutable ARRAY<int> flagged_for_removal;

    COMBINED_COLLISIONS_COUPLED_PULL_IN(RIGID_DEFORMABLE_COLLISIONS<TV>& rigid_deformable_collisions_input,const T dt_input);
    virtual ~COMBINED_COLLISIONS_COUPLED_PULL_IN();

    virtual void Discover(const T dt,const T time);
    virtual int Count() const;
    virtual void Precompute(int e);
    virtual void Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const;
    TV Compute_Impulse(int e,const COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>& ci,T scale,T dt,T time,bool diff) const;
    virtual bool Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const;
    virtual T Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const;
    T Effective_Mass(int e) const;
    virtual void Moved(int e);
    virtual T Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const { /* This isn't going to work... */ return 0; }
//#####################################################################
};
}
#endif
