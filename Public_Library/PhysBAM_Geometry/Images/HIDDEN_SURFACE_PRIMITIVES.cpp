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
SURFACE_PRIMITIVE(const TV &a,const int pa)
    :parent(pa)
{ 
    num_vertices=1;
    vertices(0)=a;
    bounding_box=RANGE<TV2>(a.Remove_Index(2));
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SURFACE_PRIMITIVE<T>::
SURFACE_PRIMITIVE(const TV &a,const TV &b,const int pa)
    :parent(pa)
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
SURFACE_PRIMITIVE(const TV &a,const TV &b,const TV &c,const int pa)
    :parent(pa)
{ 
    num_vertices=3;
    vertices(0)=a;
    vertices(1)=b;
    vertices(2)=c;
    bounding_box=RANGE<TV2>::Bounding_Box(a.Remove_Index(2),b.Remove_Index(2),c.Remove_Index(2));
}
template<class T> TRIANGLE_3D<T> SURFACE_PRIMITIVE<T>::
As_Triangle() const
{
    PHYSBAM_ASSERT(num_vertices==3);
    return TRIANGLE_3D<T>(vertices);
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
Add_Edge(HASHTABLE<VECTOR<int,2> >& edges,int a,int b)
{
    edges.Set(VECTOR<int,2>(a,b));
}
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Handle_Intersection_Triangle_Triangle(int a,int b,ARRAY<ARRAY<int> >& adjacency_list,
    ARRAY<VECTOR<int,2> >& pairs,HASHTABLE<VECTOR<int,2> >& edges)
{
    primitives.Preallocate(primitives.m+2); // Don't resize on us later

    int index[2]={a,b};
    SURFACE_PRIMITIVE<T>* primitive[2]={&primitives(index[0]),&primitives(index[1])};
    typename TRIANGLE_3D<T>::INTERSECTS_HELPER ih[2];
    TRIANGLE_3D<T> triangle[2]={primitive[0]->As_Triangle(),primitive[1]->As_Triangle()};
    bool intersects=triangle[0].Intersects(triangle[1],1e-12,ih);
    if(ih[0].n.z==0 && ih[1].n.z==0){
        T rz=(ih[0].n.Cross(ih[1].n)).z;
        if(!rz || !ih[0].neg || !ih[0].pos || !ih[1].neg || !ih[1].pos) return;
        int i=rz>0?1:-1;
        int in_front=max(i*ih[0].t[0],i*ih[0].t[1])<max(i*ih[1].t[0],i*ih[1].t[1]);
        Add_Edge(edges,index[in_front],index[1-in_front]);
        return;}

    int cut_index=abs(ih[0].n.z)>abs(ih[1].n.z);
    if(!intersects){
        bool ci_diff=ih[cut_index].neg && ih[cut_index].pos;
        bool nci_diff=ih[1-cut_index].neg && ih[1-cut_index].pos;
        if(!ci_diff || !nci_diff){
            if(ci_diff) cut_index=1-cut_index;
            int ind=ih[cut_index].w.Arg_Abs_Max();
            int in_front=((ih[cut_index].w[ind]>0)!=(ih[1-cut_index].n.z>0))^cut_index;
            Add_Edge(edges,index[in_front],index[1-in_front]);
            return;}

        int st=(ih[0].t[0]+ih[0].t[1]>ih[1].t[0]+ih[1].t[1]),ind0=argmax(ih[0].t[0],ih[0].t[1]),ind1=argmin(ih[1].t[0],ih[1].t[1]);
        if(st){ind0=1-ind0;ind1=1-ind1;}
        TV u(ih[0].x(ih[0].is)-ih[0].x(ih[0].i[ind0]),0,ih[0].w(ih[0].is)-ih[0].w(ih[0].i[ind0]));
        TV v(ih[1].x(ih[1].is)-ih[1].x(ih[1].i[ind1]),ih[1].w(ih[1].is)-ih[1].w(ih[1].i[ind1]),0);
        TV e(ih[0].n.Cross(ih[1].n).z,ih[0].n.z,ih[1].n.z);
        int in_front=(TV::Triple_Product(u,v,e)<0)^st^(u.z>0)^(v.y>0);
        Add_Edge(edges,index[in_front],index[1-in_front]);
        return;
    }

    const TV& X0=primitive[cut_index]->vertices(ih[cut_index].is);
    TV& X1=primitive[cut_index]->vertices(ih[cut_index].i[0]);
    TV& X2=primitive[cut_index]->vertices(ih[cut_index].i[1]);
    TV cut0=X0+ih[cut_index].th[0]*(X1-X0);
    TV cut1=X0+ih[cut_index].th[1]*(X2-X0);
    int new0=ih[cut_index].th[0]<1?Add(cut0,X1,X2,primitive[cut_index]->parent):-1;
    int new1=ih[cut_index].th[1]<1?Add(cut0,X2,cut1,primitive[cut_index]->parent):-1;
    X1=cut0;
    X2=cut1;

    if(new0>=0) adjacency_list.Append(ARRAY<int>());
    if(new1>=0) adjacency_list.Append(ARRAY<int>());
    ARRAY<int>& adj=adjacency_list(index[cut_index]);
    for(int i=0;i<adj.m;i++){
        if(new0>=0) adjacency_list(new0).Append(adj(i));
        if(new1>=0) adjacency_list(new1).Append(adj(i));
        if(new0>=0) adjacency_list(adj(i)).Append(new0);
        if(new1>=0) adjacency_list(adj(i)).Append(new1);
        if(adj(i)==index[1-cut_index]) continue;
        if(edges.Contains(VECTOR<int,2>(index[cut_index],adj(i)))){
            if(new0>=0) Add_Edge(edges,new0,adj(i));
            if(new1>=0) Add_Edge(edges,new1,adj(i));}
        else if(edges.Contains(VECTOR<int,2>(adj(i),index[cut_index]))){
            if(new0>=0) Add_Edge(edges,adj(i),new0);
            if(new1>=0) Add_Edge(edges,adj(i),new1);}
        else{
            if(new0>=0) pairs.Append(VECTOR<int,2>(new0,adj(i)));
            if(new1>=0) pairs.Append(VECTOR<int,2>(new1,adj(i)));}}

    int in_front=((ih[cut_index].w(ih[cut_index].is)>0)!=(ih[1-cut_index].n.z>0))^cut_index;
    Add_Edge(edges,index[in_front],index[1-in_front]);
    if(new0>=0){
        index[cut_index]=new0;
        Add_Edge(edges,index[1-in_front],index[in_front]);}
    if(new0>=1){
        index[cut_index]=new1;
        Add_Edge(edges,index[1-in_front],index[in_front]);}
}
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Handle_Intersection(int a,int b,ARRAY<ARRAY<int> >& adjacency_list,ARRAY<VECTOR<int,2> >& pairs,
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
        case 3*3: Handle_Intersection_Triangle_Triangle(a,b,adjacency_list,pairs,edges);break;}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Initialize(DIRECTED_GRAPH<>& dg)
{
    BOX_HIERARCHY<TV2> bh;
    for(int i=0;i<primitives.m;i++)
        bh.box_hierarchy.Append(primitives(i).bounding_box);
    bh.Set_Leaf_Boxes(bh.box_hierarchy,true);
    ARRAY<ARRAY<int> > intersection_list(primitives.m),adjacency_list(primitives.m);
    BOX_VISITOR_TRIVIAL visitor(intersection_list);
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
        if(Projections_Intersect(pairs(i).x,pairs(i).y))
            Handle_Intersection(pairs(i).x,pairs(i).y,adjacency_list,pairs,edges);

    dg.Initialize(primitives.m);
    for(HASHTABLE<VECTOR<int,2> >::ITERATOR it(edges);it.Valid();it.Next())
        if(Test_Edge(it.Key().x,it.Key().y))
            dg.Add_Edge(it.Key().x,it.Key().y);
}
template<class T> bool HIDDEN_SURFACE_PRIMITIVES<T>::
Divide_Primitive(int divide,int cutter,ARRAY<int>& inside,ARRAY<int>& outside)
{
    TRIANGLE_3D<T> triangle(primitives(divide).As_Triangle());
    PLANE<T> plane0=Get_Cutting_Plane(primitives(cutter).vertices(0),primitives(cutter).vertices(1));
    PLANE<T> plane1=Get_Cutting_Plane(primitives(cutter).vertices(1),primitives(cutter).vertices(2));
    PLANE<T> plane2=Get_Cutting_Plane(primitives(cutter).vertices(2),primitives(cutter).vertices(0));
    bool flip=!plane0.Lazy_Inside(primitives(cutter).vertices(2));
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
    for(int i=1;i<inside_triangles.m;i++)
        inside.Append(Add(inside_triangles(i).x1,inside_triangles(i).x2,inside_triangles(i).x3,primitives(divide).parent));
    for(int i=0;i<outside_triangles.m;i++)
        outside.Append(Add(outside_triangles(i).x1,outside_triangles(i).x2,outside_triangles(i).x3,primitives(divide).parent));

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
template<class T> int HIDDEN_SURFACE_PRIMITIVES<T>::
Add(const TV &a,int pa)
{
    return primitives.Append(SURFACE_PRIMITIVE<T>(a,pa>=0?pa:primitives.m));
}
template<class T> int HIDDEN_SURFACE_PRIMITIVES<T>::
Add(const TV &a,const TV &b,int pa)
{
    return primitives.Append(SURFACE_PRIMITIVE<T>(a,b,pa>=0?pa:primitives.m));
}
template<class T> int HIDDEN_SURFACE_PRIMITIVES<T>::
Add(const TV &a,const TV &b,const TV &c,int pa)
{
    return primitives.Append(SURFACE_PRIMITIVE<T>(a,b,c,pa>=0?pa:primitives.m));
}
template class HIDDEN_SURFACE_PRIMITIVES<float>;
template class HIDDEN_SURFACE_PRIMITIVES<double>;
