//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE
//#####################################################################
#ifndef __COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE__
#define __COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE__
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>

namespace PhysBAM{

template<class TV> class TRIANGLE_COLLISIONS;
template<class TV> class TRIANGLE_REPULSIONS;
template<class TV> struct COMBINED_COLLISIONS_DEFORMABLE_IMPULSE;

template<class TV,int d>
class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE:public COMBINED_COLLISIONS<TV>::COLLIDER
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename COMBINED_COLLISIONS<TV>::IMPULSE IMPULSE;
public:

    struct COLLISION
    {
        VECTOR<T,d> weights;
        TV normal;
        VECTOR<int,d> nodes;
        TV relative_velocity;
        T target_relative_velocity;
        mutable TV total_impulse;
    };

    TRIANGLE_COLLISIONS<TV>& triangle_collisions;
    TRIANGLE_REPULSIONS<TV>& triangle_repulsions;
    ARRAY<COLLISION> collisions;
    bool update_swept_hierarchies;
    mutable ARRAY<bool> flagged_for_removal;

    COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE(TRIANGLE_COLLISIONS<TV>& triangle_collisions_input,TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,const bool update_swept_hierarchies_input);
    virtual ~COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE();

    virtual int Count() const;
    virtual void Precompute(int e);
    virtual void Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const;
    virtual bool Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const;
    virtual T Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const;
    virtual void Moved(int e);
    virtual T Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const;
    T Effective_Mass(int e) const;
    TV Compute_Impulse(int e,const COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>& di,T scale,T dt,T time,bool diff) const;
//#####################################################################
};
}
#endif
