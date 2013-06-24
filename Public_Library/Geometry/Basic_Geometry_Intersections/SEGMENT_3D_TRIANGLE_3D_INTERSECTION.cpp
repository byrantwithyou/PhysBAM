//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const SEGMENT_3D<T>& segment,const TRIANGLE_3D<T>& triangle,const T thickness_over_two)
{
    RAY<VECTOR<T,3> > ray(segment.X.x,segment.X.y-segment.X.x);
    ray.semi_infinite=false;
    ray.t_max=(segment.X.y-segment.X.x).Magnitude();
    return INTERSECTION::Intersects(ray,triangle,thickness_over_two);
}
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const SEGMENT_3D<T>& segment,const TRIANGLE_3D<T>& triangle,VECTOR<T,2>& weights_segment,VECTOR<T,3>& weights,const T thickness_over_two)
{
    RAY<VECTOR<T,3> > ray(segment.X.x,segment.X.y-segment.X.x);
    ray.semi_infinite=false;
    T magnitude=(segment.X.y-segment.X.x).Magnitude();
    ray.t_max=magnitude;
    if(INTERSECTION::Intersects(ray,triangle,thickness_over_two)){
        weights_segment(0)=ray.t_max/magnitude;
        weights_segment(1)=1-weights_segment(0);
        weights=triangle.Barycentric_Coordinates(ray.Point(ray.t_max));
        return true;}
    return false;
}
//#####################################################################
template bool Intersects(const SEGMENT_3D<float>&,const TRIANGLE_3D<float>&,const float);
template bool Intersects(const SEGMENT_3D<float>&,const TRIANGLE_3D<float>&,VECTOR<float,2>&,VECTOR<float,3>&,const float);
template bool Intersects(const SEGMENT_3D<double>&,const TRIANGLE_3D<double>&,const double);
template bool Intersects(const SEGMENT_3D<double>&,const TRIANGLE_3D<double>&,VECTOR<double,2>&,VECTOR<double,3>&,const double);
};
};
