//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SEGMENT_ORIGIN_AREAS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_ORIGIN_AREAS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/VOLUME_COLLISIONS.h>
#include <algorithm>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
template<class TV>
struct VOLUME_COLLISIONS_VISITOR
{
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_OBJECT;
    ARRAY<VECTOR<int,2> > pairs;
    const T_OBJECT& ta1;
    const T_OBJECT& ta2;

    VOLUME_COLLISIONS_VISITOR(const T_OBJECT& a,const T_OBJECT& b): ta1(a),ta2(b) {}

    ~VOLUME_COLLISIONS_VISITOR(){}

    bool Cull_Self(const int a) const
    {return false;}

    bool Cull(const int a,const int b) const
    {return false;}

    void Store(const int a,const int b)
    {if(Topology_Aware_Intersection_Test(ta1.mesh.elements(a),ta2.mesh.elements(b),(ARRAY_VIEW<const TV>)ta1.particles.X)) pairs.Append(VECTOR<int,2>(a,b));}
};
//#####################################################################
// Function Compute_Collision_Triangles
//#####################################################################
template<class TV> void VOLUME_COLLISIONS<TV>::
Compute_Collision_Triangles()
{
    area=0;
    gradient.Remove_All();
    hessian.Remove_All();
    for(int i=0;i<objects.m;i++) objects(i)->Initialize_Hierarchy();
    for(int i=0;i<objects.m;i++)
        for(int j=i;j<=objects.m;j++)
            Compute_Collision_Triangles(*objects(i),*objects(j));
}
int compute_collision_triangles_order[2][4][3]=
{
    {{2,3},{3,1},{1,2}},
    {{1,2,3},{2,3,4},{3,4,1},{4,1,2}}
};
static int Sort_Pair(VECTOR<int,4>& I)
{
    int sign=1;
    if(I(1)<I(2)){exchange(I(1),I(2));sign=-sign;}
    if(I(3)<I(4)){exchange(I(3),I(4));sign=-sign;}
    if(I(1)<I(3)) for(int k=0;k<2;k++) exchange(I(k),I(k+2));
    return sign;
}
static int Sort_Pair(VECTOR<int,6>& I)
{
    int sign=1;
    if(I(1)<I(2)){exchange(I(1),I(2));sign=-sign;}
    if(I(1)<I(3)){exchange(I(1),I(3));sign=-sign;}
    if(I(2)<I(3)){exchange(I(2),I(3));sign=-sign;}

    if(I(4)<I(5)){exchange(I(4),I(5));sign=-sign;}
    if(I(4)<I(6)){exchange(I(4),I(6));sign=-sign;}
    if(I(5)<I(6)){exchange(I(5),I(6));sign=-sign;}

    if(I(1)<I(4)) for(int k=0;k<3;k++) exchange(I(k),I(k+3));
    return sign;
}
//#####################################################################
// Function Compute_Collision_Triangles
//#####################################################################
template<class TV> void VOLUME_COLLISIONS<TV>::
Compute_Collision_Triangles(T_OBJECT& obj1,T_OBJECT& obj2)
{
    VOLUME_COLLISIONS_VISITOR<TV> visitor(obj1,obj2);
    obj1.hierarchy->Intersection_List(*obj2.hierarchy,visitor,ZERO());
    if(visitor.pairs.m==0) return;

    HASHTABLE<int> visited_particles1;
    HASHTABLE<int> visited_particles2;
    HASHTABLE<VECTOR<int,TV::m*2>,int> to_process;
    ARRAY_VIEW<TV> X(obj1.particles.X);

    for(int m=0;m<visitor.pairs.m;m++){
        int a,b;visitor.pairs(m).Get(a,b);
        //if(TRIANGLE_2D<T>::Signed_Size(X.Subset(obj1.mesh.elements(a)))
        //  *TRIANGLE_2D<T>::Signed_Size(X.Subset(obj2.mesh.elements(b)))<=0)
        //    continue;
        visited_particles1.Set_All(obj1.mesh.elements(a));
        visited_particles2.Set_All(obj2.mesh.elements(b));
        for(int i=0;i<TV::m+1;i++){
            for(int j=0;j<TV::m+1;j++){
                VECTOR<int,TV::m*2> I;
                for(int k=0;k<TV::m;k++){
                    I(k)=obj1.mesh.elements(a)(compute_collision_triangles_order[TV::m-2][i-1][k-1]);
                    I(k+TV::m)=obj2.mesh.elements(b)(compute_collision_triangles_order[TV::m-2][j-1][k-1]);}
                int sign=Sort_Pair(I);
                to_process.Get_Or_Insert(I)+=sign;}}}

    TV x0;
    for(HASHTABLE<int>::ITERATOR it(visited_particles1);it.Valid();it.Next())
        x0+=X(it.Key());
    for(HASHTABLE<int>::ITERATOR it(visited_particles2);it.Valid();it.Next())
        x0+=X(it.Key());
    x0/=static_cast<T>(visited_particles1.Size()+visited_particles2.Size());

    for(typename HASHTABLE<VECTOR<int,TV::m*2>,int>::ITERATOR it(to_process);it.Valid();it.Next()){
        T sign=static_cast<T>(it.Data());
        if(sign==0) continue;
        VECTOR<int,TV::m*2> I=it.Key();
        ORIGIN_AREAS::VOL_DATA<T,TV::m,TV::m*2> data;
        TV PTS[TV::m*2];
        for(int k=1;k<=TV::m*2;k++) PTS[k-1]=X(I(k));
        Volume_From_Simplices(data,x0,PTS); // gotta love ADL!
        area+=sign*data.V;
        for(int k=1;k<=TV::m*2;k++) gradient.Get_Or_Insert(I(k))+=sign*data.G[k-1];
        for(int k=1;k<=TV::m*2;k++) for(int m=1;m<=TV::m*2;m++) hessian.Get_Or_Insert(VECTOR<int,2>(I(k),I(m)))+=sign*data.H[k-1][m-1];}
}
template class VOLUME_COLLISIONS<VECTOR<float,2> >;
template class VOLUME_COLLISIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VOLUME_COLLISIONS<VECTOR<double,2> >;
template class VOLUME_COLLISIONS<VECTOR<double,3> >;
#endif
