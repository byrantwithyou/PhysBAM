//#####################################################################
// Copyright 2005-2007, Kevin Der, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOLUME_COLLISIONS
//##################################################################### 
#ifndef __VOLUME_COLLISIONS__
#define __VOLUME_COLLISIONS__    

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
namespace PhysBAM{

template<class TV>
class VOLUME_COLLISIONS
{
    typedef typename TV::SCALAR T;
public:
    ARRAY<TRIANGULATED_AREA<T>*> triangulated_areas;

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

    void Compute_Collision_Edges();
    void Compute_Collision_Edges(TRIANGULATED_AREA<T>& ta,SEGMENTED_CURVE_2D<T>& sc);
};   
}
#endif

