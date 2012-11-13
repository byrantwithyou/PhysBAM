//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_2D_SEGMENT_2D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(const SEGMENT_2D<T>& segment1,const SEGMENT_2D<T>& segment2,const T thickness_over_two)
{
    RAY<VECTOR<T,2> > ray(segment2.X.x,segment2.X.y-segment2.X.x);ray.semi_infinite=false;ray.t_max=(segment2.X.y-segment2.X.x).Magnitude();
    return INTERSECTION::Intersects(ray,segment1,thickness_over_two);
}
//#####################################################################
template bool Intersects(const SEGMENT_2D<float>&,const SEGMENT_2D<float>&,const float);
template bool Intersects(const SEGMENT_2D<double>&,const SEGMENT_2D<double>&,const double);
};
};
