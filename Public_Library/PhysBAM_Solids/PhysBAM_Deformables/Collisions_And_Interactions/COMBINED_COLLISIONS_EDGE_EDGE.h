//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_EDGE_EDGE
//#####################################################################
#ifndef __COMBINED_COLLISIONS_EDGE_EDGE__
#define __COMBINED_COLLISIONS_EDGE_EDGE__
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_REPULSION_BASE.h>

namespace PhysBAM{

template<class TV> class TRIANGLE_REPULSIONS;

template<class TV>
class COMBINED_COLLISIONS_EDGE_EDGE:public COMBINED_COLLISIONS_REPULSION_BASE<TV>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename COMBINED_COLLISIONS<TV>::IMPULSE IMPULSE;
    typedef COMBINED_COLLISIONS_REPULSION_BASE<TV> BASE;
    using BASE::collisions;using BASE::prune_pairs;using BASE::triangle_repulsions;
public:

    COMBINED_COLLISIONS_EDGE_EDGE(TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,bool prune_pairs_input);
    virtual ~COMBINED_COLLISIONS_EDGE_EDGE();

    virtual void Discover(const T dt,const T time);
    void Discover_Pruned();
    void Discover_Saved_Pairs();
//#####################################################################
};
}
#endif
