//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_POINT_FACE_PULL_IN
//#####################################################################
#ifndef __COMBINED_COLLISIONS_POINT_FACE_PULL_IN__
#define __COMBINED_COLLISIONS_POINT_FACE_PULL_IN__
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE.h>

namespace PhysBAM{

template<class TV> class TRIANGLE_COLLISIONS;
template<class TV> class TRIANGLE_REPULSIONS;

template<class TV>
class COMBINED_COLLISIONS_POINT_FACE_PULL_IN:public COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,TV::m+1>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
    typedef typename COMBINED_COLLISIONS<TV>::IMPULSE IMPULSE;
    typedef COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d+1> BASE;
    using BASE::collisions;using BASE::COLLISION;using BASE::triangle_collisions;using BASE::triangle_repulsions;
    using BASE::update_swept_hierarchies;using BASE::flagged_for_removal;

    HASHTABLE<VECTOR<int,d+1>,bool> point_face_pairs;
public:
    const T dt;


    COMBINED_COLLISIONS_POINT_FACE_PULL_IN(TRIANGLE_COLLISIONS<TV>& triangle_collisions_input,TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,const T dt_input,const bool update_swept_hierarchies_input);
    virtual ~COMBINED_COLLISIONS_POINT_FACE_PULL_IN();

    virtual void Discover(const T dt,const T time);
    void Discover_Pruned();
    //void Discover_Saved_Pairs();
//#####################################################################
};
}
#endif
