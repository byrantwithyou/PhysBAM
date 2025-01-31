//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __RAY_TETRAHEDRON_INTERSECTION__
#define __RAY_TETRAHEDRON_INTERSECTION__
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const TETRAHEDRON<T>& tetrahedron,const T thickness=0);
//#####################################################################
};
};
#endif
