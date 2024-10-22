//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Intersections/RAY_ORIENTED_BOX_INTERSECTION.h>
#include <Geometry/Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Fuzzy_Intersects(RAY<VECTOR<T,2> >& ray,const ORIENTED_BOX<VECTOR<T,2> >& box,const T segment_intersect_epsilon)
{
    bool intersects=false;
    for(int i=0;i<2;i++) if(INTERSECTION::Fuzzy_Intersects(ray,SEGMENT_2D<T>(box.corner,box.corner+box.edges.Column(i)),segment_intersect_epsilon)) intersects=true;
    for(int i=0;i<2;i++) if(INTERSECTION::Fuzzy_Intersects(ray,SEGMENT_2D<T>(box.corner+box.edges.Column(i),box.corner+box.edges.Column_Sum()),segment_intersect_epsilon)) intersects=true;
    return intersects;
}
//#####################################################################
template bool Fuzzy_Intersects(RAY<VECTOR<float,2> >&,const ORIENTED_BOX<VECTOR<float,2> >&,const float);
template bool Fuzzy_Intersects(RAY<VECTOR<double,2> >&,const ORIENTED_BOX<VECTOR<double,2> >&,const double);
};
};
