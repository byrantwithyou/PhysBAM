//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_INTERSECTION
//##################################################################### 
#ifndef __TRIANGLE_INTERSECTION__
#define __TRIANGLE_INTERSECTION__    

#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
namespace PhysBAM{
template<class T,class TV> T Triangle_Intersection_Area(const TRIANGLE_2D<T>& a,const TRIANGLE_2D<T>& b,VECTOR<TV,6>& G,VECTOR<VECTOR<MATRIX<T,2>,6>,6>& H);
template<class TV> bool Topology_Aware_Triangle_Intersection_Test(VECTOR<int,3> a,VECTOR<int,3> b,ARRAY_VIEW<const TV> X);
}
#endif

