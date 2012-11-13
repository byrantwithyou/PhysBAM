//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TETRAHEDRON_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
// Function Intersects
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const TETRAHEDRON<T>& tetrahedron, const T thickness)
{
    bool intersection=false;
    for(int i=0;i<4;i++)
        if(INTERSECTION::Intersects(ray,tetrahedron.triangle(i),thickness)){
            intersection=true;
            ray.aggregate_id=i;}
    return intersection;
}
//#####################################################################
template bool Intersects(RAY<VECTOR<float,3> >&,const TETRAHEDRON<float>&,const float);
template bool Intersects(RAY<VECTOR<double,3> >&,const TETRAHEDRON<double>&,const double);
};
};
