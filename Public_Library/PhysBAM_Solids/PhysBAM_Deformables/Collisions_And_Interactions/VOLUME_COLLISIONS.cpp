//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SEGMENT_ORIGIN_AREAS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/VOLUME_COLLISIONS.h>
#include <algorithm>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
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
    for(int i=1;i<=triangulated_areas.m;i++) triangulated_areas(i)->Initialize_Hierarchy();
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
    VOLUME_COLLISIONS_VISITOR<T> visitor(ta1,ta2);
    ta1.hierarchy->Intersection_List(*ta2.hierarchy,visitor,ZERO());

    HASHTABLE<VECTOR<int,4>,int> to_process;
    ARRAY_VIEW<TV> X(ta1.particles.X);

    for(int i=1;i<=visitor.pairs.m;i++){
        int a,b;visitor.pairs(i).Get(a,b);

        VECTOR<int,6> me;
        ta1.mesh.elements(a).Get(me(1),me(2),me(3));
        ta2.mesh.elements(b).Get(me(4),me(5),me(6));

        for(int i=1;i<=3;i++){
            int in=i%3+1;
            for(int j=1;j<=3;j++){
                int jn=j%3+1;
                VECTOR<int,4> I(me(i),me(in),me(j+3),me(jn+3));
                int sign=1;
                if(I(1)<I(2)){exchange(I(1),I(2));sign=-sign;}
                if(I(3)<I(4)){exchange(I(3),I(4));sign=-sign;}
                if(I(1)<I(3)){exchange(I(1),I(3));exchange(I(2),I(4));}
                to_process.Get_Or_Insert(I)+=sign;}}}

    for(typename HASHTABLE<VECTOR<int,4>,int>::ITERATOR it(to_process);it.Valid();it.Next()){
        T sign=it.Data();
        if(sign==0) continue;
        VECTOR<int,4> I=it.Key();
        SEGMENT_ORIGIN_AREAS::DATA<T,1,4> data;
        Area_From_Segments(data,X(I(1)),X(I(2)),X(I(3)),X(I(4)));
        area+=sign*data.V.x;
        for(int k=1;k<=4;k++) gradient.Get_Or_Insert(I(k))+=sign*(const TV&)data.G[k-1];
        for(int k=1;k<=4;k++) for(int m=1;m<=4;m++) hessian.Get_Or_Insert(VECTOR<int,2>(I(k),I(m)))+=sign*data.H[0][k-1][m-1];}
}
template class VOLUME_COLLISIONS<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VOLUME_COLLISIONS<VECTOR<double,2> >;
#endif
