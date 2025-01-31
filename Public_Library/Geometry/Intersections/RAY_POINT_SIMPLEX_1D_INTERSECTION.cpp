//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <Core/Vectors/VECTOR_1D.h>
#include <Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Intersections/RAY_BOX_INTERSECTION.h>
#include <Geometry/Intersections/RAY_POINT_SIMPLEX_1D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,1> >& ray,const POINT_SIMPLEX_1D<T>& point,const T thickness_over_two)
{
    return Intersects(ray,RANGE<VECTOR<T,1> >(point.X.x),thickness_over_two);
}
//#####################################################################
// Function Closest_Non_Intersecting_Point
//#####################################################################
template<class T> bool Closest_Non_Intersecting_Point(RAY<VECTOR<T,1> >& ray,const POINT_SIMPLEX_1D<T>& point,const T thickness_over_two)
{
    RANGE<VECTOR<T,1> > thickened_box=RANGE<VECTOR<T,1> >::Bounding_Box(point.X.x-2*thickness_over_two,point.X.x+2*thickness_over_two);
    return Intersects(ray,RANGE<VECTOR<T,1> >(point.X.x),thickness_over_two);
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,1> >&,const POINT_SIMPLEX_1D<float>&,const float);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<float,1> >&,const POINT_SIMPLEX_1D<float>&,const float);
template bool Intersects(RAY<VECTOR<double,1> >&,const POINT_SIMPLEX_1D<double>&,const double);
template bool Closest_Non_Intersecting_Point(RAY<VECTOR<double,1> >&,const POINT_SIMPLEX_1D<double>&,const double);
};
};
