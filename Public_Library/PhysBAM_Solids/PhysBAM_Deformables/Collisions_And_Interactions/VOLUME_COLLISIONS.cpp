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
    const SEGMENTED_CURVE<TV>& segments;

    VOLUME_COLLISIONS_VISITOR(const SEGMENTED_CURVE<TV>& segments_input): segments(segments_input) {}

    ~INTERSECTING_PAIRS_VISITOR(){}

    bool Cull_Self(const int a) const
    {return false;}

    bool Cull(const int a,const int b) const
    {return false;}

    void Store(const int a,const int b)
    {VECTOR<int,2> A(segments.mesh.elements(a)),B(segments.mesh.elements(b));
    if(A.x!=B.x && A.x!=B.y && A.y!=B.x && A.y!=B.y) candidate_pairs.Append(VECTOR<int,2>(a,b));}
//#####################################################################
};
}
//#####################################################################
// Function Initialize_Components
//#####################################################################
template<class TV> void VOLUME_COLLISIONS<TV>::
Initialize_Components()
{
    segments.mesh.Initialize_Adjacent_Elements();
    component_lookup.Resize(segments.mesh.elements.m);

    for(int i=1;i<=segments.mesh.elements.m;i++)
        if(!component_lookup(i)){
            int ci=components.Append(COMPONENT());
            COMPONENT& c=components(ci);
            c.closed=true;
            int first=i,prev=first;
            component_lookup(first)=ci;
            c.list.Append(first);
            const ARRAY<int>& first_adj=segments.mesh.adjacent_elements(first);
            PHYSBAM_ASSERT(first_adj.m && first_adj<=2);
            for(int cur=first_adj(1);cur!=first;){
                component_lookup(cur)=ci;
                c.list.Append(cur);
                const ARRAY<int>& adj=segments.mesh.adjacent_elements(cur);
                PHYSBAM_ASSERT(adj.m && adj<=2);
                if(adj.m==1){
                    if(!c.closed) break;
                    c.closed=false;
                    if(first_adj.m==1) break;
                    std::reverse(c.list.begin(),c.list.end());
                    prev=first;
                    cur=first_adj(2);}
                else{
                    int p=prev;
                    prev=cur;
                    cur=adj(1)^adj(2)^p;}}
}
//#####################################################################
// Function Compute_Loops
//#####################################################################
template<class TV> void VOLUME_COLLISIONS<TV>::
Compute_Loops()
{
    if(!segments.hierarchy) segments.Initialize_Hierarchy();

    VOLUME_COLLISIONS_VISITOR<TV> visitor;
    segments.hierarchy.Intersection_List(segments.hierarchy,visitor);

    HASHTABLE<VECTOR<int,2>,ARRAY<VECTOR<int,2> > > component_intersections;

    for(int i=1;i<=visitor.m;i++){
        int a,b;visitor(i).Get(a,b);
        if(Intersects(segments.Get_Element(a),segments.Get_Element(b),0)){
            component_intersections.Get_Or_Insert(VECTOR<int,2>(component_lookup(a),component_lookup(b))).Append(VECTOR<int,2>(a,b));}}

    for(HASHTABLE<VECTOR<int,2>,ARRAY<int,2> >::ITERATOR it(component_intersections);it.Valid();it.Next()){
        PHYSBAM_ASSERT(it->Key().x==it->Key().y || it->Data().m%2==0);}

    
}

