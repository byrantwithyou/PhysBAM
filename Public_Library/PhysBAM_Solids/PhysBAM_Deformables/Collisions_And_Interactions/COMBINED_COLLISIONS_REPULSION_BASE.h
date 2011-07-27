//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_REPULSION_BASE
//#####################################################################
#ifndef __COMBINED_COLLISIONS_REPULSION_BASE__
#define __COMBINED_COLLISIONS_REPULSION_BASE__
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>

namespace PhysBAM{

template<class TV> class TRIANGLE_REPULSIONS;

template<class TV>
class COMBINED_COLLISIONS_REPULSION_BASE:public COMBINED_COLLISIONS<TV>::COLLIDER
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename COMBINED_COLLISIONS<TV>::IMPULSE IMPULSE;
public:

    struct COLLISION
    {
        VECTOR<T,4> weights;
        TV normal;
        VECTOR<int,4> nodes;
        TV relative_velocity;
    };

    TRIANGLE_REPULSIONS<TV>& triangle_repulsions;
    ARRAY<COLLISION> collisions;

    bool prune_pairs;

    COMBINED_COLLISIONS_REPULSION_BASE(TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,bool prune_pairs_input);
    virtual ~COMBINED_COLLISIONS_REPULSION_BASE();

    virtual int Count() const;
    virtual void Precompute(int e);
    virtual void Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const;
    virtual bool Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const;
    virtual T Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const;
    virtual void Moved(int e);
    virtual T Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const;
    T Effective_Mass(int e) const;
    TV Compute_Impulse(int e,T scale,T dt,T time,bool diff) const;
//#####################################################################
};
}
#endif
