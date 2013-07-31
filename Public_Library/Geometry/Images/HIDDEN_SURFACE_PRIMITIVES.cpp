#ifdef USE_BOOST_GEOMETRY
//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Tools/Log/LOG.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Images/HIDDEN_SURFACE_PRIMITIVES.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <boost/geometry/algorithms/centroid.hpp>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> SURFACE_PRIMITIVE<T>::
SURFACE_PRIMITIVE(const TRIANGLE_3D<T> &t,const int pa)
    :triangle(t),parent(pa)
{ 
}
template<class T> void SURFACE_PRIMITIVE<T>::
Init_Projection()
{
    projection.Resize(1);
    projection(0).outer().Append(triangle.X.x.Remove_Index(2));
    projection(0).outer().Append(triangle.X.y.Remove_Index(2));
    projection(0).outer().Append(triangle.X.z.Remove_Index(2));
    projection(0).outer().Append(triangle.X.x.Remove_Index(2));
    boost::geometry::correct(projection);
}
template<class T> bool HIDDEN_SURFACE_PRIMITIVES<T>::
Projections_Intersect(int a,int b)
{
    typename SURFACE_PRIMITIVE<T>::MULTI_POLYGON in;
    boost::geometry::intersection(primitives(a).projection,primitives(b).projection,in);
    return boost::geometry::area(in)>1e-10;
//    return boost::geometry::intersects(primitives(a).projection,primitives(b).projection);
}
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Add_Edge(HASHTABLE<VECTOR<int,2> >& edges,int a,int b)
{
    edges.Set(VECTOR<int,2>(a,b));
}
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Add_Edge(HASHTABLE<VECTOR<int,2> >& edges,int a,int b,TV2 pt)
{
    TV Na=primitives(a).triangle.Normal();
    TV Nb=primitives(b).triangle.Normal();
    T za=-(pt.Append(0)-primitives(a).triangle.X.x).Dot(Na)/Na.z;
    T zb=-(pt.Append(0)-primitives(b).triangle.X.x).Dot(Nb)/Nb.z;
    if(za>zb) Add_Edge(edges,a,b);
    else Add_Edge(edges,b,a);
}
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Handle_Intersection(int a,int b,ARRAY<ARRAY<int> >& adjacency_list,ARRAY<VECTOR<int,2> >& pairs,
    HASHTABLE<VECTOR<int,2> >& edges)
{
    T area_tol=1e-10;
    primitives.Preallocate(primitives.m+1); // Don't resize on us later

    typename SURFACE_PRIMITIVE<T>::MULTI_POLYGON projected_intersection;
    boost::geometry::intersection(primitives(a).projection,primitives(b).projection,projected_intersection);
    if(boost::geometry::area(projected_intersection)<area_tol) return;

    int index[2]={a,b};
    SURFACE_PRIMITIVE<T>* primitive[2]={&primitives(index[0]),&primitives(index[1])};
    typename TRIANGLE_3D<T>::INTERSECTS_HELPER ih[2];
    TRIANGLE_3D<T> triangle[2]={primitive[0]->triangle,primitive[1]->triangle};
    bool intersects=triangle[0].Intersects(triangle[1],1e-12,ih);
    if(!intersects){
        TV2 pt;
        boost::geometry::centroid(projected_intersection,pt);
        return Add_Edge(edges,a,b,pt);}

    TV2 X0=primitive[0]->triangle.X(ih[0].is).Remove_Index(2);
    TV2 X1=primitive[0]->triangle.X(ih[0].i[0]).Remove_Index(2);
    TV2 X2=primitive[0]->triangle.X(ih[0].i[1]).Remove_Index(2);
    TV2 cut0=X0+ih[0].th[0]*(X1-X0);
    TV2 cut1=X0+ih[0].th[1]*(X2-X0);

    typename SURFACE_PRIMITIVE<T>::MULTI_POLYGON tri,quad,int_tri,int_quad;
    tri.Resize(1);
    tri(0).outer().Append(X0);
    tri(0).outer().Append(cut0);
    tri(0).outer().Append(cut1);
    tri(0).outer().Append(X0);
    boost::geometry::correct(tri);
    quad.Resize(1);
    quad(0).outer().Append(cut0);
    quad(0).outer().Append(X1);
    quad(0).outer().Append(X2);
    quad(0).outer().Append(cut1);
    quad(0).outer().Append(cut0);
    boost::geometry::correct(quad);

    boost::geometry::intersection(tri,primitives(a).projection,int_tri);
    boost::geometry::intersection(quad,primitives(a).projection,int_quad);
    T area_tri=boost::geometry::area(int_tri);
    T area_quad=boost::geometry::area(int_quad);
    if(area_tri<area_tol && area_quad<area_tol)
        return;

    if(area_tri<area_tol){
        int_quad.Exchange(primitives(a).projection);
        TV2 pt;
        boost::geometry::centroid(primitives(a).projection,pt);
        return Add_Edge(edges,a,b,pt);}

    int_tri.Exchange(primitives(a).projection);
    TV2 pt;
    boost::geometry::centroid(primitives(a).projection,pt);
    Add_Edge(edges,a,b,pt);
    if(area_quad<area_tol) return;

    int outside=Add(primitives(a).triangle,primitives(a).parent);
    int_quad.Exchange(primitives(outside).projection);
    boost::geometry::centroid(primitives(outside).projection,pt);
    Add_Edge(edges,outside,b,pt);

    adjacency_list.Append(ARRAY<int>());
    ARRAY<int>& adj=adjacency_list(a);
    for(int i=0;i<adj.m;i++){
        adjacency_list(outside).Append(adj(i));
        adjacency_list(adj(i)).Append(outside);
        if(adj(i)==b) continue;
        if(edges.Contains(VECTOR<int,2>(a,adj(i)))) Add_Edge(edges,outside,adj(i));
        else if(edges.Contains(VECTOR<int,2>(adj(i),a))) Add_Edge(edges,adj(i),outside);
        else pairs.Append(VECTOR<int,2>(outside,adj(i)));}
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void HIDDEN_SURFACE_PRIMITIVES<T>::
Initialize(DIRECTED_GRAPH<>& dg)
{
    BOX_HIERARCHY<TV2> bh;
    for(int i=0;i<primitives.m;i++){
        bh.box_hierarchy.Append(primitives(i).triangle.Bounding_Box().Remove_Dimension(2));
        primitives(i).Init_Projection();}
    bh.Set_Leaf_Boxes(bh.box_hierarchy,true);
    ARRAY<ARRAY<int> > intersection_list(primitives.m),adjacency_list(primitives.m);
    BOX_VISITOR_TRIVIAL visitor(intersection_list);
    bh.Intersection_List(bh,visitor,0);

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
        Handle_Intersection(pairs(i).x,pairs(i).y,adjacency_list,pairs,edges);

    dg.Initialize(primitives.m);
    for(HASHTABLE<VECTOR<int,2> >::ITERATOR it(edges);it.Valid();it.Next())
        if(Test_Edge(it.Key().x,it.Key().y))
            dg.Add_Edge(it.Key().x,it.Key().y);
}
template<class T> int HIDDEN_SURFACE_PRIMITIVES<T>::
Divide_Primitive(int divide,int cutter)
{
    typename SURFACE_PRIMITIVE<T>::MULTI_POLYGON in,out;
    boost::geometry::intersection(primitives(divide).projection,primitives(cutter).projection,in);
    boost::geometry::difference(primitives(divide).projection,primitives(cutter).projection,out);
    if(!out.Size()) return -1;
    int outside=Add(primitives(divide).triangle,primitives(divide).parent);
    out.Exchange(primitives(outside).projection);
    return outside;
}
template<class T> int HIDDEN_SURFACE_PRIMITIVES<T>::
Add(const TRIANGLE_3D<T> &t,int pa)
{
    return primitives.Append(SURFACE_PRIMITIVE<T>(t,pa>=0?pa:primitives.m));
}
namespace PhysBAM{
template class HIDDEN_SURFACE_PRIMITIVES<float>;
template class HIDDEN_SURFACE_PRIMITIVES<double>;
}
#endif
