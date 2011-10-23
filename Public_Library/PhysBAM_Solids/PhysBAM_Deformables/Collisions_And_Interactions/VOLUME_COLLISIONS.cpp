//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/VOLUME_COLLISIONS.h>
#include <algorithm>
using namespace PhysBAM;

template<class T>
struct VOLUME_COLLISIONS_VISITOR
{
    ARRAY<VECTOR<int,2> > pairs;
    const TRIANGULATED_AREA<T>& ta1;
    const TRIANGULATED_AREA<T>& ta2;

    VOLUME_COLLISIONS_VISITOR(const TRIANGULATED_AREA<T>& a,const TRIANGULATED_AREA<T>& b): ta1(a),ta2(b) {}

    ~VOLUME_COLLISIONS_VISITOR(){}

    bool Cull_Self(const int a) const
    {return false;}

    bool Cull(const int a,const int b) const
    {return false;}

    void Store(const int a,const int b)
    {if(Topology_Aware_Triangle_Intersection_Test(ta1.mesh.elements(a),ta2.mesh.elements(b),(ARRAY_VIEW<const VECTOR<T,2> >)ta1.particles.X)) pairs.Append(VECTOR<int,2>(a,b));}

//#####################################################################
};
//#####################################################################
// Function Compute_Loops
//#####################################################################
template<class T> void VOLUME_COLLISIONS<VECTOR<T,2> >::
Compute_Collision_Triangles()
{
    area=0;
    gradient.Remove_All();
    hessian.Remove_All();
    for(int i=1;i<=triangulated_areas.m;i++)
        for(int j=i;j<=triangulated_areas.m;j++)
            Compute_Collision_Triangles(*triangulated_areas(i),*triangulated_areas(j));
}
//#####################################################################
// Function Compute_Loops
//#####################################################################
template<class T> void VOLUME_COLLISIONS<VECTOR<T,2> >::
Compute_Collision_Triangles(TRIANGULATED_AREA<T>& ta1,TRIANGULATED_AREA<T>& ta2)
{
    if(!ta1.hierarchy) ta1.Initialize_Hierarchy();
    if(!ta2.hierarchy) ta2.Initialize_Hierarchy();

    VOLUME_COLLISIONS_VISITOR<T> visitor(ta1,ta2);
    ta1.hierarchy->Intersection_List(*ta2.hierarchy,visitor,ZERO());

    for(int i=1;i<=visitor.pairs.m;i++){
        int a,b;visitor.pairs(i).Get(a,b);
        VECTOR<TV,6> G;
        VECTOR<VECTOR<MATRIX<T,2>,6>,6> H;
        area+=Triangle_Intersection_Area(ta1.Get_Element(a),ta2.Get_Element(b),G,H);
        VECTOR<int,6> I;
        ta1.mesh.elements(a).Get(I(1),I(2),I(3));
        ta2.mesh.elements(b).Get(I(4),I(5),I(6));
        for(int k=1;k<=6;k++) gradient.Get_Or_Insert(I(k))=G(k);
        for(int k=1;k<=6;k++) for(int m=1;m<=6;m++) hessian.Get_Or_Insert(VECTOR<int,2>(I(k),I(m)))+=H(k)(m);}
}
template class VOLUME_COLLISIONS<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VOLUME_COLLISIONS<VECTOR<double,2> >;
#endif
