//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_RIGID
//#####################################################################
#ifndef __COMBINED_COLLISIONS_RIGID__
#define __COMBINED_COLLISIONS_RIGID__
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>

namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLISIONS;
template<class TV> class RIGID_BODY_COLLECTION;

// TODO: Set last checked
// TODO: flags: ignore separating, collision thickness, etc.
// TODO: pairs_processed_by_collisions

template<class TV>
class COMBINED_COLLISIONS_RIGID:public COMBINED_COLLISIONS<TV>::COLLIDER
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename COMBINED_COLLISIONS<TV>::IMPULSE IMPULSE;
public:

    struct COLLISION
    {
        int particle_body;
        int particle_index;
        int levelset_body;
        TV location;
        TV normal;
        TV relative_velocity;
    };

    RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions;
    ARRAY<COLLISION> collisions;

    bool ignore_separating;
    bool rolling_friction;
    bool clamp_energy;
    bool clamp_friction_magnitude;
    T desired_separation_distance;
    bool use_coefficient_of_restitution;
    bool use_pairs_processed_by_collisions;
    bool use_saved_pairs;
    bool use_collisions;

    COMBINED_COLLISIONS_RIGID(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions_input,bool do_collision,bool use_saved_pairs_input);
    virtual ~COMBINED_COLLISIONS_RIGID();

    virtual void Discover(const T dt,const T time);
    virtual int Count() const;
    virtual void Precompute(int e);
    virtual void Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const;
    virtual bool Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const;
    virtual T Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const;
    virtual void Moved(int e);
    virtual T Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const;

    void Contact_Discover();
    void Collision_Discover(const T dt,const T time);
    void Discover_Helper(const ARRAY<VECTOR<int,2> >& pairs);
    TWIST<TV> Compute_Impulse_Derivative(int e,T scale,const T coefficient_of_restitution,const T coefficient_of_friction,const bool clamp_friction_magnitude) const;
//#####################################################################
};
}
#endif
