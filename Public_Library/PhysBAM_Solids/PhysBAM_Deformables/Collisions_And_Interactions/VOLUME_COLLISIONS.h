//#####################################################################
// Copyright 2005-2007, Kevin Der, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOLUME_COLLISIONS
//##################################################################### 
#ifndef __VOLUME_COLLISIONS__
#define __VOLUME_COLLISIONS__    

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV>
class VOLUME_COLLISIONS
{
public:
    SEGMENTED_CURVE<TV> segments;

    struct COMPONENT
    {
        ARRAY<int> list;
        bool closed;
    };

    ARRAY<COMPONENT> components;
    ARRAY<int> component_lookup;

    struct LOOP
    {
        VECTOR<ARRAY<int>,2> parts;
    };

    ARRAY<LOOP> loops;

    void Initialize_Components();
    void Compute_Loops();

};   
}
#endif

