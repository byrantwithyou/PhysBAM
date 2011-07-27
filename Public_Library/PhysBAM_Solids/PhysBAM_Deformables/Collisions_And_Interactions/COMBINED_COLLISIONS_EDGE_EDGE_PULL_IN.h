//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN
//#####################################################################
#ifndef __COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN__
#define __COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE.h>

namespace PhysBAM{

template<class TV> class TRIANGLE_COLLISIONS;
template<class TV> class TRIANGLE_REPULSIONS;

template<class TV>
class COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN:public COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,2*TV::m-2>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
    typedef typename COMBINED_COLLISIONS<TV>::IMPULSE IMPULSE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,d-1>::SIMPLEX_FACE T_EDGE;
    typedef COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,2*d-2> BASE;
    using BASE::collisions;using BASE::triangle_collisions;using BASE::triangle_repulsions;
    using BASE::update_swept_hierarchies;using BASE::flagged_for_removal;

    HASHTABLE<VECTOR<int,2*d-2>,bool> edge_edge_pairs;
public:
    const T dt;

    COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN(TRIANGLE_COLLISIONS<TV>& triangle_collisions_input,TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,const T dt,const bool update_swept_hierarchies_input);
    virtual ~COMBINED_COLLISIONS_EDGE_EDGE_PULL_IN();

    virtual void Discover(const T dt,const T time);
    void Discover_Pruned();
    void Discover_Saved_Pairs();
//#####################################################################
};
}
#endif
