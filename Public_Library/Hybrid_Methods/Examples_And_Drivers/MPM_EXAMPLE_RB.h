//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_EXAMPLE_RB__
#define __MPM_EXAMPLE_RB__
#include <Core/Data_Structures/CHAINED_ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_COLLECTION;

template<class TV>
class MPM_EXAMPLE_RB:public MPM_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_EXAMPLE_RB(const STREAM_TYPE stream_type_input);
    MPM_EXAMPLE_RB(const MPM_EXAMPLE_RB&) = delete;
    void operator=(const MPM_EXAMPLE_RB&) = delete;
    virtual ~MPM_EXAMPLE_RB();

    // First entry MUST be int, and its value MUST be nonnegative.
    struct RASTERIZED_DATA
    {
        int id;
        T phi;
    };
    CHAINED_ARRAY<RASTERIZED_DATA,TV_INT> rasterized_data;
    ARRAY<bool> rigid_body_is_simulated;
    HASHTABLE<PAIR<int,int>,TV> stored_contacts_rr;
    HASHTABLE<PAIR<int,TV_INT>,TV> stored_contacts_rm;
    T contact_factor=1; // omega
    T impulse_interpolation=1; // lambda
    T min_impulse_change=1e-4;
    
    bool pairwise_collisions=false;
    bool projected_collisions=false;
    int collision_iterations=5;
//#####################################################################
};
}
#endif
