//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/ELLIPSOID.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Intersections/RAY_ELLIPSOID_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const ELLIPSOID<VECTOR<T,3> >& ellipsoid,const T thickness)
{ // TODO(jontg): everything else takes thickness_over_two; this should be made consistent...
    typedef VECTOR<T,3> TV;
    // Apply inverse transformation to ray
    TV ray_endpoint=ellipsoid.radii.Inverse_Times(ellipsoid.orientation.Inverse_Rotate(ray.endpoint-ellipsoid.center));
    TV ray_direction=ellipsoid.radii.Inverse_Times(ellipsoid.orientation.Inverse_Rotate(ray.direction)).Normalized();
    
    // sphere intersection routine
    T thickness_over_two=(T).5*thickness;
    T distance_squared=TV::Dot_Product(ray_endpoint,ray_endpoint);
    T outside_shell_squared=sqr((T)1+thickness_over_two);
    T c=distance_squared-outside_shell_squared;
    if(c > 0){ // outside   
        T b=TV::Dot_Product(ray_direction,ray_endpoint);if(b >= 0) return false; // no intersection - ray goes away from sphere 
        T d=sqr(b)-c;if(d < 0) return 0; // no intersection - ray misses sphere
        T t=-b-sqrt(d); // smaller of the two roots
        if(ray.semi_infinite || t < ray.t_max){ray.semi_infinite=false;ray.t_max=t;return true;}
        else return false;}
    else{ // inside or on the boundary
        T inside_shell_squared=sqr((T)1-thickness_over_two);
        T c=distance_squared-inside_shell_squared;
        if(c < 0){ // inside
            T b=TV::Dot_Product(ray_direction,ray_endpoint);
            T d=sqr(b)-c;
            T t=-b+sqrt(d); // larger of the two roots
            if(ray.semi_infinite || t < ray.t_max){ray.semi_infinite=false;ray.t_max=t;return true;}
            else return false;}
        else{ray.semi_infinite=false;ray.t_max=0;return true;}} // on the boundary
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const ELLIPSOID<VECTOR<float,3> >&,const float);
template bool Intersects(RAY<VECTOR<double,3> >&,const ELLIPSOID<VECTOR<double,3> >&,const double);
};
};
