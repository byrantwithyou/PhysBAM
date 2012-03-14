//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Images/HIDDEN_SURFACE_PRIMITIVES.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> SURFACE_PRIMITIVE<T>::
SURFACE_PRIMITIVE(const TV &a,const int pa):parent(pa)
{ 
    num_vertices=1;
    vertices(0)=a;
    bounding_box=RANGE<TV2>::Bounding_Box(a.Remove_Index(2));
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SURFACE_PRIMITIVE<T>::
SURFACE_PRIMITIVE(const TV &a,const TV &b,const int pa):parent(pa)
{ 
    num_vertices=2;
    vertices(0)=a;
    vertices(1)=b;
    bounding_box=RANGE<TV2>::Bounding_Box(a.Remove_Index(2),b.Remove_Index(2));
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SURFACE_PRIMITIVE<T>::
SURFACE_PRIMITIVE(const TV &a,const TV &b,const TV &c,const int pa):parent(pa)
{ 
    num_vertices=3;
    vertices(0)=a;
    vertices(1)=b;
    vertices(2)=c;
    bounding_box=RANGE<TV2>::Bounding_Box(a.Remove_Index(2),b.Remove_Index(2),c.Remove_Index(2));
}
template<class T> bool HIDDEN_SURFACE_PRIMITIVES<T>::
Projections_Intersect(int a,int b)
{
    const SURFACE_PRIMITIVE<T>& pa=primitives(a),&pb=primitives(b);
    if(pa.num_vertices>pb.num_vertices) return Projections_Intersect(b,a);
    switch(pa.num_vertices*pb.num_vertices){
        case 1*1: PHYSBAM_FATAL_ERROR();break;
        case 1*2: PHYSBAM_FATAL_ERROR();break;
        case 1*3: PHYSBAM_FATAL_ERROR();break;
        case 2*2: PHYSBAM_FATAL_ERROR();break;
        case 2*3: PHYSBAM_FATAL_ERROR();break;
        case 3*3:
            TRIANGLE_2D<T> ta(pa.vertices(0).Remove_Index(2),pa.vertices(1).Remove_Index(2),pa.vertices(2).Remove_Index(2));
            TRIANGLE_2D<T> tb(pb.vertices(0).Remove_Index(2),pb.vertices(1).Remove_Index(2),pb.vertices(2).Remove_Index(2));
            return ta.Intersects(tb);
    }
    PHYSBAM_FATAL_ERROR();
}
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Add_Edge(DIRECTED_GRAPH<>& dg,HASHTABLE<VECTOR<int,2> >& edges,int a,int b)
{
    dg.Add_Edge(a,b);
    edges.Set(VECTOR<int,2>(a,b));
}
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Handle_Intersection_Triangle_Triangle(DIRECTED_GRAPH<>& dg,int a,int b,ARRAY<ARRAY<int> >& adjacency_list,
    ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& edges)
{
    int index[2]={a,b};
    SURFACE_PRIMITIVE<T>* primitive[2]={&primitives(index[0]),&primitives(index[1])};
    typename TRIANGLE_3D<T>::INTERSECTS_HELPER ih[2];
    TRIANGLE_3D<T> triangle[2]={TRIANGLE_3D<T>(primitive[0]->vertices),TRIANGLE_3D<T>(primitive[1]->vertices)};
    bool intersects=triangle[0].Intersects(triangle[1],ih);
    if(ih[0].n.z==0 && ih[1].n.z==0){
        T rz=(ih[0].n.Cross(ih[1].n)).z;
        if(!rz || !ih[0].neg || !ih[0].pos || !ih[1].neg || !ih[1].pos) return;
        int i=rz>0?1:-1;
        int in_front=max(i*ih[0].t[0],i*ih[0].t[1])<max(i*ih[1].t[0],i*ih[1].t[1]);
        Add_Edge(dg,edges,index[in_front],index[1-in_front]);
        return;}

    int cut_index=abs(ih[0].n.z)>abs(ih[1].n.z);
    if(!intersects){
        int ind=ih[cut_index].w.Arg_Abs_Max();
        int in_front=(ih[cut_index].w[ind]>0)==(ih[1-cut_index].n.z>0);
        Add_Edge(dg,edges,index[in_front],index[1-in_front]);
        return;}

    const TV& X0=primitive[cut_index]->vertices(ih[cut_index].i[0]);
    TV& X1=primitive[cut_index]->vertices(ih[cut_index].i[1]);
    TV& X2=primitive[cut_index]->vertices(ih[cut_index].i[2]);
    TV cut0=X0+ih[cut_index].t[0]*(X1-X0);
    TV cut1=X0+ih[cut_index].t[1]*(X2-X0);
    int new0=primitives.Append(SURFACE_PRIMITIVE<T>(cut0,X1,X2,primitive[cut_index]->parent));
    int new1=primitives.Append(SURFACE_PRIMITIVE<T>(cut0,X2,cut1,primitive[cut_index]->parent));
    X1=cut0;
    X2=cut1;

    int adj0=adjacency_list.Append(ARRAY<int>());
    int adj1=adjacency_list.Append(ARRAY<int>());
    PHYSBAM_ASSERT(adj0==new0 && adj1==new1);
    ARRAY<int>& adj=adjacency_list(index[cut_index]);
    for(int i=0;i<adj.m;i++){
        adjacency_list(new0).Append(adj(i));
        adjacency_list(new1).Append(adj(i));
        adjacency_list(adj(i)).Append(new0);
        adjacency_list(adj(i)).Append(new1);
        if(adj(i)==index[1-cut_index]) continue;
        if(edges.Contains(VECTOR<int,2>(index[cut_index],adj(i)))){
            Add_Edge(dg,edges,new0,adj(i));
            Add_Edge(dg,edges,new1,adj(i));}
        else if(edges.Contains(VECTOR<int,2>(adj(i),index[cut_index]))){
            Add_Edge(dg,edges,adj(i),new0);
            Add_Edge(dg,edges,adj(i),new1);}
        else{
            pairs.Append(VECTOR<int,2>(new0,adj(i)));
            pairs.Append(VECTOR<int,2>(new1,adj(i)));}}

    int in_front=(ih[cut_index].w[ih[cut_index].i[0]]>0)==(ih[1-cut_index].n.z>0);
    Add_Edge(dg,edges,index[in_front],index[1-in_front]);
    index[cut_index]=new0;
    Add_Edge(dg,edges,index[1-in_front],index[in_front]);
    index[cut_index]=new1;
    Add_Edge(dg,edges,index[1-in_front],index[in_front]);
}
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Handle_Intersection(DIRECTED_GRAPH<>& dg,int a,int b,ARRAY<ARRAY<int> >& adjacency_list,ARRAY<VECTOR<int,2> >& pairs,
    HASHTABLE<VECTOR<int,2> >& edges)
{
    int index[2]={a,b};
    SURFACE_PRIMITIVE<T>* primitive[2]={&primitives(index[0]),&primitives(index[1])};
    switch(primitive[0]->num_vertices*primitive[1]->num_vertices){
        case 1*1: PHYSBAM_FATAL_ERROR();break;
        case 1*2: PHYSBAM_FATAL_ERROR();break;
        case 1*3: PHYSBAM_FATAL_ERROR();break;
        case 2*2: PHYSBAM_FATAL_ERROR();break;
        case 2*3: PHYSBAM_FATAL_ERROR();break;
        case 3*3: Handle_Intersection_Triangle_Triangle(dg,a,b,adjacency_list,pairs,edges);break;}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Initialize(DIRECTED_GRAPH<>& dg)
{
    dg.Initialize(primitives.m);
    BOX_HIERARCHY<TV2> bh;
    for(int i=0;i<primitives.m;i++)
        bh.box_hierarchy.Append(primitives(i).bounding_box);
    bh.Set_Leaf_Boxes(bh.box_hierarchy,true);
    ARRAY<ARRAY<int> > intersection_list(primitives.m),adjacency_list(primitives.m);
    BOX_VISITOR_TRIVIAL visitor(intersection_list);
//    bh.Intersection_List(visitor);

    bh.Intersection_List(bh,visitor,ZERO());

    HASHTABLE<VECTOR<int,2> > edges;
    ARRAY<VECTOR<int,2> > pairs;
    for(int i=0;i<intersection_list.m;i++)
        for(int j=0;j<intersection_list(i).m;j++){
            int k=intersection_list(i)(j);
            if(Projections_Intersect(i,k)){
                pairs.Append(VECTOR<int,2>(i,k));
                adjacency_list(i).Append(k);
                adjacency_list(k).Append(i);}}

    for(int i=0;i<pairs.m;i++)
        Handle_Intersection(dg,pairs(i).x,pairs(i).y,adjacency_list,pairs,edges);
}
template<class T> bool HIDDEN_SURFACE_PRIMITIVES<T>::
Divide_Primitive(int divide,int cutter,ARRAY<int>& inside,ARRAY<int>& outside)
{
    TRIANGLE_3D<T> triangle(primitives(divide).vertices);
    PLANE<T> plane0=Get_Cutting_Plane(primitives(cutter).vertices(0),primitives(cutter).vertices(1));
    PLANE<T> plane1=Get_Cutting_Plane(primitives(cutter).vertices(1),primitives(cutter).vertices(2));
    PLANE<T> plane2=Get_Cutting_Plane(primitives(cutter).vertices(2),primitives(cutter).vertices(0));
    bool flip=plane0.Lazy_Inside(primitives(cutter).vertices(2));
    if(flip){
        plane0.normal=-plane0.normal;
        plane1.normal=-plane1.normal;
        plane2.normal=-plane2.normal;}
    ARRAY<TRIANGLE_3D<T> > inside_triangles0,inside_triangles1,inside_triangles,outside_triangles;
    triangle.Cut_With_Hyperplane(triangle,plane0,inside_triangles0,outside_triangles);
    for(int i=0;i<inside_triangles0.m;i++)
        triangle.Cut_With_Hyperplane(inside_triangles0(i),plane1,inside_triangles1,outside_triangles);
    for(int i=0;i<inside_triangles1.m;i++)
        triangle.Cut_With_Hyperplane(inside_triangles1(i),plane2,inside_triangles,outside_triangles);

    if(!inside_triangles.m) return false; // Cutting wasn't necessary
    primitives(divide).vertices(0)=inside_triangles(0).x1;
    primitives(divide).vertices(1)=inside_triangles(0).x2;
    primitives(divide).vertices(2)=inside_triangles(0).x3;
    inside.Append(divide);
    for(int i=1;i<inside_triangles.m;i++){
        SURFACE_PRIMITIVE<T> sp(inside_triangles(i).x1,inside_triangles(i).x2,inside_triangles(i).x3,primitives(divide).parent);
        inside.Append(primitives.Append(sp));}
    for(int i=0;i<outside_triangles.m;i++){
        SURFACE_PRIMITIVE<T> sp(outside_triangles(i).x1,outside_triangles(i).x2,outside_triangles(i).x3,primitives(divide).parent);
        outside.Append(primitives.Append(sp));}

    return true;
}
//#####################################################################
// Function Get_Cutting_Plane
//#####################################################################
template<class T> PLANE<T> HIDDEN_SURFACE_PRIMITIVES<T>::
Get_Cutting_Plane(const TV &a,const TV &b)
{ 
    return PLANE<T>((a-b).Remove_Index(2).Rotate_Clockwise_90().Normalized().Append(0),a);
}
template class HIDDEN_SURFACE_PRIMITIVES<float>;
template class HIDDEN_SURFACE_PRIMITIVES<double>;
