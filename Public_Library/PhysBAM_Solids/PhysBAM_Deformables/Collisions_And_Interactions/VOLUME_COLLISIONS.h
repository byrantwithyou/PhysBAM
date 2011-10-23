//#####################################################################
// Copyright 2005-2007, Kevin Der, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOLUME_COLLISIONS
//##################################################################### 
#ifndef __VOLUME_COLLISIONS__
#define __VOLUME_COLLISIONS__    

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
namespace PhysBAM{

template<class TV> class VOLUME_COLLISIONS;

template<class T>
class VOLUME_COLLISIONS<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
public:
    ARRAY<TRIANGULATED_AREA<T>*> triangulated_areas;
    T area;
    HASHTABLE<int,TV> gradient;
    HASHTABLE<VECTOR<int,2>,MATRIX<T,2> > hessian;

    void Compute_Collision_Triangles();
    void Compute_Collision_Triangles(TRIANGULATED_AREA<T>& ta1,TRIANGULATED_AREA<T>& ta2);
};   
}
#endif

