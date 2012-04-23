//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __RAY_TRIANGLE_3D_INTERSECTION__
#define __RAY_TRIANGLE_3D_INTERSECTION__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T,class TV> bool Intersects(RAY<TV>& ray,const TRIANGLE_3D<T>& plane,const T thickness_over_two=0);
template<class T,class TV> bool Lazy_Intersects(RAY<TV>& ray,const TRIANGLE_3D<T>& plane,const T thickness_over_two=0);
template<class T,class TV> bool Closest_Non_Intersecting_Point(RAY<TV>& ray,const TRIANGLE_3D<T>& plane,const T thickness_over_two=0);
//#####################################################################
};
};
#endif
