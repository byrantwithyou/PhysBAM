//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Intersections/RAY_PLANE_INTERSECTION.h>
#include <Geometry/Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T,class TV> bool Intersects(RAY<TV>& ray,const TRIANGLE_3D<T>& triangle,const T thickness_over_two)
{
    RAY<TV> ray_temp;
    ray.Save_Intersection_Information(ray_temp);
    T thickness=2*thickness_over_two;

    // first check the plane of the triangle's face
    TV normal=triangle.Normal();
    if(!INTERSECTION::Intersects(ray,PLANE<T>(normal,triangle.X.x),thickness_over_two))
        return false; // otherwise intersects the plane
    TV plane_point(ray.Point(ray.t_max));
    PLANE<T> edge_plane_12(TV::Cross_Product(triangle.X.y-triangle.X.x,normal).Normalized(),triangle.X.x),edge_plane_23(TV::Cross_Product(triangle.X.z-triangle.X.y,normal).Normalized(),triangle.X.y),
                     edge_plane_31(TV::Cross_Product(triangle.X.x-triangle.X.z,normal).Normalized(),triangle.X.z);
    if(!edge_plane_12.Outside(plane_point,thickness) && !edge_plane_23.Outside(plane_point,thickness) && !edge_plane_31.Outside(plane_point,thickness))return true; // intersects face of triangle wedge 
    else ray.Restore_Intersection_Information(ray_temp);

    PLANE<T> triangle_plane(normal,triangle.X.x);
    // check for intersection with the sides of the wedge
    if(edge_plane_12.Outside(plane_point,thickness_over_two) && INTERSECTION::Intersects(ray,edge_plane_12,thickness_over_two)){
        TV edge_point(ray.Point(ray.t_max));
        if(triangle_plane.Boundary(edge_point,thickness) && !edge_plane_23.Outside(edge_point,thickness) && !edge_plane_31.Outside(edge_point,thickness)){
            ray.intersection_location=RAY<TV>::INTERIOR_POINT;return true;}
        else ray.Restore_Intersection_Information(ray_temp);}
    if(edge_plane_23.Outside(plane_point,thickness_over_two) && INTERSECTION::Intersects(ray,edge_plane_23,thickness_over_two)){
        TV edge_point(ray.Point(ray.t_max));
        if(triangle_plane.Boundary(edge_point,thickness) && !edge_plane_12.Outside(edge_point,thickness) && !edge_plane_31.Outside(edge_point,thickness)){
            ray.intersection_location=RAY<TV>::INTERIOR_POINT;return true;}
        else ray.Restore_Intersection_Information(ray_temp);}
    if(edge_plane_31.Outside(plane_point,thickness_over_two) && INTERSECTION::Intersects(ray,edge_plane_31,thickness_over_two)){
        TV edge_point(ray.Point(ray.t_max));
        if(triangle_plane.Boundary(edge_point,thickness) && !edge_plane_12.Outside(edge_point,thickness) && !edge_plane_23.Outside(edge_point,thickness)){
            ray.intersection_location=RAY<TV>::INTERIOR_POINT;return true;}
        else ray.Restore_Intersection_Information(ray_temp);}

    return false;
}
//#####################################################################
// Function Lazy_Intersects
//#####################################################################
template<class T,class TV> bool Lazy_Intersects(RAY<TV>& ray,const TRIANGLE_3D<T>& triangle,const T thickness_over_two)
{
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate_id=ray.aggregate_id;
    PLANE<T> triangle_plane(triangle.Normal(),triangle.X.x);
    if(!INTERSECTION::Intersects(ray,triangle_plane)) return false; // otherwise intersects the plane
    if(triangle.Lazy_Planar_Point_Inside_Triangle(ray.Point(ray.t_max))) return true; // intersects the face of the triangle 
    else{ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate_id;return false;} // reset ray
}
//#####################################################################
// Function Closest_Non_Intersecting_Point
//#####################################################################
template<class T,class TV> bool Closest_Non_Intersecting_Point(RAY<TV>& ray,const TRIANGLE_3D<T>& triangle,const T thickness_over_two)
{
    RAY<TV> ray_temp;
    ray.Save_Intersection_Information(ray_temp);
    if(!INTERSECTION::Intersects(ray,triangle,thickness_over_two)) return false;
    if(ray.intersection_location==RAY<TV>::START_POINT) return true;
    ray.Restore_Intersection_Information(ray_temp);

    TV normal=triangle.Normal();
    // TODO: Save having to re-generate all the planes...
    T thickness=2*thickness_over_two;
    TV normal_times_thickness=normal*thickness;
    TRIANGLE_3D<T> top_triangle(triangle.X+normal_times_thickness);
    TRIANGLE_3D<T> bottom_triangle(triangle.X-normal_times_thickness);
    PLANE<T> top_plane(normal,triangle.X.x);
    PLANE<T> bottom_plane(normal,triangle.X.x);
    PLANE<T> edge_plane_12(TV::Cross_Product(triangle.X.y-triangle.X.x,normal).Normalized(),triangle.X.x);
    edge_plane_12.x0+=edge_plane_12.normal*thickness;
    PLANE<T> edge_plane_23(TV::Cross_Product(triangle.X.z-triangle.X.y,normal).Normalized(),triangle.X.y);
    edge_plane_23.x0+=edge_plane_23.normal*thickness;
    PLANE<T> edge_plane_31(TV::Cross_Product(triangle.X.x-triangle.X.z,normal).Normalized(),triangle.X.z);
    edge_plane_31.x0+=edge_plane_31.normal*thickness;
    bool found_intersection=false;
    if(INTERSECTION::Intersects(ray,top_triangle,thickness_over_two)) found_intersection=true;
    if(INTERSECTION::Intersects(ray,bottom_triangle,thickness_over_two)) found_intersection=true;
    if(INTERSECTION::Rectangle_Intersects(ray,edge_plane_12,top_plane,bottom_plane,edge_plane_23,edge_plane_31,thickness_over_two)) found_intersection=true;
    if(INTERSECTION::Rectangle_Intersects(ray,edge_plane_23,top_plane,bottom_plane,edge_plane_12,edge_plane_31,thickness_over_two)) found_intersection=true;
    if(INTERSECTION::Rectangle_Intersects(ray,edge_plane_31,top_plane,bottom_plane,edge_plane_12,edge_plane_23,thickness_over_two)) found_intersection=true;
    return found_intersection;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const TRIANGLE_3D<float>&,const float);
template bool Lazy_Intersects(RAY<VECTOR<float,3> >&,const TRIANGLE_3D<float>&,const float);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<float,3> >&,const TRIANGLE_3D<float>&, const float);
template bool Intersects(RAY<VECTOR<double,3> >&,const TRIANGLE_3D<double>&,const double);
template bool Lazy_Intersects(RAY<VECTOR<double,3> >&,const TRIANGLE_3D<double>&,const double);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<double,3> >&,const TRIANGLE_3D<double>&, const double);
};
};
