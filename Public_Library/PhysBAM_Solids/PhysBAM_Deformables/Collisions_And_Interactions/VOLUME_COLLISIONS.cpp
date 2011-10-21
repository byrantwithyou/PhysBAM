//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/VOLUME_COLLISIONS.h>
#include <algorithm>
using namespace PhysBAM;

namespace{
template<class TV>
struct VOLUME_COLLISIONS_VISITOR
{
    typedef typename TV::SCALAR T;
    ARRAY<VECTOR<int,2> > candidate_pairs;
    const TRIANGULATED_AREA<T>& ta;
    const SEGMENTED_CURVE<T>& sc;

//    VOLUME_COLLISIONS_VISITOR(const SEGMENTED_CURVE<TV>& segments_input): segments(segments_input) {}

    ~VOLUME_COLLISIONS_VISITOR(){}

    bool Cull_Self(const int a) const
    {return false;}

    bool Cull(const int s,const int t) const
    {return false;}

    void Store(const int s,const int t)
    {if(!ta.mesh.elements(t).Contains_All(sc.mesh.elements(s))) candidate_pairs.Append(VECTOR<int,2>(s,t));}
//#####################################################################
};
}
//#####################################################################
// Function Compute_Loops
//#####################################################################
template<class TV> void VOLUME_COLLISIONS<TV>::
Compute_Collision_Edges()
{
    for(int i=1;i<=triangulated_areas;i++)
        for(int j=1;j<=triangulated_areas;j++)
            Compute_Collision_Edges(*triangulated_areas(i),triangulated_areas(i)->Get_Boundary_Object());
}
//#####################################################################
// Function Compute_Loops
//#####################################################################
template<class TV> void VOLUME_COLLISIONS<TV>::
Compute_Collision_Edges(TRIANGULATED_AREA<T>& ta,SEGMENTED_CURVE_2D<T>& sc)
{
    if(!ta.hierarchy) ta.Initialize_Hierarchy();
    if(!sc.hierarchy) sc.Initialize_Hierarchy();

    VOLUME_COLLISIONS_VISITOR<TV> visitor;
    sc.hierarchy.Intersection_List(ta.hierarchy,visitor);
/*
    for(int i=1;i<=visitor.m;i++){
        int a,b;visitor(i).Get(a,b);
        if(Intersects(sc.Get_Element(a),ta.Get_Element(b),0)){
            component_intersections.Get_Or_Insert(VECTOR<int,2>(component_lookup(a),component_lookup(b))).Append(VECTOR<int,2>(a,b));}}

    for(HASHTABLE<VECTOR<int,2>,ARRAY<int,2> >::ITERATOR it(component_intersections);it.Valid();it.Next()){
        PHYSBAM_ASSERT(it->Key().x==it->Key().y || it->Data().m%2==0);}
*/
    
}



