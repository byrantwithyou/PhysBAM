//#####################################################################
// Copyright 2005-2007, Kevin Der, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOLUME_COLLISIONS
//##################################################################### 
#ifndef __VOLUME_COLLISIONS__
#define __VOLUME_COLLISIONS__    

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV>
class VOLUME_COLLISIONS
{
    typedef typename TV::SCALAR T;typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_OBJECT;
public:
    ARRAY<T_OBJECT*> objects;
    T area;
    HASHTABLE<int,TV> gradient;
    HASHTABLE<VECTOR<int,2>,MATRIX<T,TV::m> > hessian;

    void Compute_Collision_Triangles();
    void Compute_Collision_Triangles(T_OBJECT& obj1,T_OBJECT& obj2);
};   
}
#endif

