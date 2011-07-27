//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_COUPLED
//#####################################################################
#ifndef __COMBINED_COLLISIONS_COUPLED__
#define __COMBINED_COLLISIONS_COUPLED__
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_DEFORMABLE_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/COMBINED_COLLISIONS_RIGID.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/COMBINED_COLLISIONS_COUPLED_IMPULSE.h>

namespace PhysBAM{

template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class RIGID_DEFORMABLE_COLLISIONS;

// TODO: Update spatial partition.

template<class TV>
class COMBINED_COLLISIONS_COUPLED:public virtual COMBINED_COLLISIONS<TV>::COLLIDER
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename COMBINED_COLLISIONS<TV>::IMPULSE IMPULSE;
public:
    RIGID_DEFORMABLE_COLLISIONS<TV>& rigid_deformable_collisions;

    struct COLLISION
    {
        int body;
        int p;
        T depth;
        TV normal;
        TV relative_velocity;
        COLLISION_GEOMETRY_ID collision_body_id;
    };

    ARRAY<COLLISION> collisions;
    T collision_body_thickness;
    bool use_particles_collided_with_rigid_body;
    bool use_particles_contacting_rigid_body;
    bool clamp_friction_magnitude;
    bool use_collisions;
    bool use_saved_pairs;
    ARRAY_VIEW<const TV> X_save;

    COMBINED_COLLISIONS_COUPLED(RIGID_DEFORMABLE_COLLISIONS<TV>& rigid_deformable_collisions_input,ARRAY_VIEW<const TV> X_save_input,bool do_collision,bool use_saved_pairs_input);
    virtual ~COMBINED_COLLISIONS_COUPLED();

    virtual void Discover(const T dt,const T time);
    virtual int Count() const;
    virtual void Precompute(int e);
    virtual void Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const;
    virtual bool Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const;
    virtual T Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const;
    virtual void Moved(int e);
    virtual T Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const;
    void Contact_Discover();
    void Collision_Discover();
    void Discover_Helper(int b,int p,ARRAY_VIEW<const TV> X,COLLISION_GEOMETRY_ID collision_geometry_id);
    TV Compute_Impulse_Derivative(int e,T scale,const T coefficient_of_restitution,const T coefficient_of_friction,const bool clamp_friction_magnitude) const;
//#####################################################################
};
}
#endif
