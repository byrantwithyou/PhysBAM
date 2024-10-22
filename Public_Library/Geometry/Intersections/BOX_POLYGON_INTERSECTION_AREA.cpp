//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <Core/Data_Structures/PAIR.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/exchange.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Geometry/Basic_Geometry/POLYGON.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Intersections/BOX_POLYGON_INTERSECTION_AREA.h>
#include <Geometry/Intersections/RAY_BOX_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersection_Area
//#####################################################################
namespace {
template<class T> T Intersection_Area_Helper(const RANGE<VECTOR<T,1> >& box, const POLYGON<VECTOR<T,1> >& polygon)
{
    RANGE<VECTOR<T,1> > polygon_as_range(polygon.X(0),polygon.X(1));
    if(polygon_as_range.min_corner.x >= polygon_as_range.max_corner.x) exchange(polygon_as_range.min_corner,polygon_as_range.max_corner);
    return box.Intersection_Area(polygon_as_range);
}

template<class T> T Intersection_Area_Helper(const RANGE<VECTOR<T,2> >& box, const POLYGON<VECTOR<T,2> >& polygon)
{
    POLYGON<VECTOR<T,2> > box_polygon(box);
    POLYGON<VECTOR<T,2> > modified_polygon(polygon.X.Size());
    for(int i=0;i<polygon.X.Size();i++) modified_polygon.X(i) = polygon.X(i);

    int offset=0;int side;
    T t_start,t_end;
    ARRAY<PAIR<int,int> > added_nodes;
    for(int s=0;s<polygon.X.Size()-1;++s){
        RAY<VECTOR<T,2> > edge_ray(SEGMENT_2D<T>(polygon.X(s),polygon.X(s+1)));
        if(INTERSECTION::Get_Intersection_Range(edge_ray,box,t_start,t_end)){
            if(t_start != 0){modified_polygon.X.Insert(box_polygon.Find_Closest_Point_On_Polygon(edge_ray.Point(t_start),side),s+offset);added_nodes.Append(PAIR<int,int>(s+offset++,side));}
            if(t_end != edge_ray.t_max){modified_polygon.X.Insert(box_polygon.Find_Closest_Point_On_Polygon(edge_ray.Point(t_end),side),s+offset);added_nodes.Append(PAIR<int,int>(s+offset++,side));}}}

    RAY<VECTOR<T,2> > edge_ray(SEGMENT_2D<T>(polygon.X(polygon.X.Size()),polygon.X(0))); // Final edge, to close the polygon
    if(INTERSECTION::Get_Intersection_Range(edge_ray,box,t_start,t_end)){
        if(t_start != 0){modified_polygon.X.Append(box_polygon.Find_Closest_Point_On_Polygon(edge_ray.Point(t_start),side));added_nodes.Append(PAIR<int,int>(modified_polygon.X.Size(),side));offset++;}
        if(t_end != edge_ray.t_max){modified_polygon.X.Append(box_polygon.Find_Closest_Point_On_Polygon(edge_ray.Point(t_end),side));added_nodes.Append(PAIR<int,int>(modified_polygon.X.Size(),side));offset++;}}

    if(offset==0) { // entirely outside, or entirely inside
        if(box.Inside(polygon.X(0),0)) return modified_polygon.Area();
        else if(modified_polygon.Inside_Polygon(box.min_corner)) return box.Robust_Size();
        else return (T)0;}

    POLYGON<VECTOR<T,2> > new_polygon;
    for(int curr_intersection_node=0;curr_intersection_node<added_nodes.Size();curr_intersection_node++){
        int next_intersection_node = (curr_intersection_node % added_nodes.Size()) + 1;
        int curr_node=added_nodes(curr_intersection_node).x;
        new_polygon.X.Append(modified_polygon.X(curr_node));
        int next_node=(curr_node % modified_polygon.X.Size()) + 1;
        bool on_polygon=box.Inside(modified_polygon.X(next_node),0);
        if(on_polygon)
            for(; next_node != added_nodes(next_intersection_node).x; next_node = (next_node % modified_polygon.X.Size()) + 1)
                new_polygon.X.Append(modified_polygon.X(next_node));
        else
            for(int next_side = added_nodes(curr_intersection_node).y;next_side != added_nodes(next_intersection_node).y;next_side = (next_side % box_polygon.X.Size()) + 1)
                new_polygon.X.Append(box_polygon.X(next_side+1));}

    return new_polygon.Area();
}

template<class T> T Intersection_Area_Helper(const RANGE<VECTOR<T,3> >& box, const POLYGON<VECTOR<T,3> >& polygon)
{PHYSBAM_NOT_IMPLEMENTED();}
};
template<class TV> typename TV::SCALAR Intersection_Area(const RANGE<TV>& box, const POLYGON<TV>& polygon)
{
    return Intersection_Area_Helper<typename TV::SCALAR>(box,polygon);
}
//#####################################################################
template float Intersection_Area(const RANGE<VECTOR<float,1> >&,const POLYGON<VECTOR<float,1> >&);
template float Intersection_Area(const RANGE<VECTOR<float,2> >&,const POLYGON<VECTOR<float,2> >&);
template float Intersection_Area(const RANGE<VECTOR<float,3> >&,const POLYGON<VECTOR<float,3> >&);
template double Intersection_Area(const RANGE<VECTOR<double,1> >&,const POLYGON<VECTOR<double,1> >&);
template double Intersection_Area(const RANGE<VECTOR<double,2> >&,const POLYGON<VECTOR<double,2> >&);
template double Intersection_Area(const RANGE<VECTOR<double,3> >&,const POLYGON<VECTOR<double,3> >&);
};
};
