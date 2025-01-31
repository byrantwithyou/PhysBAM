//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __RAY_ORIENTED_BOX_INTERSECTION__
#define __RAY_ORIENTED_BOX_INTERSECTION__
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Lazy_Intersects(RAY<VECTOR<T,2> >& ray,const ORIENTED_BOX<VECTOR<T,2> >& box,const T segment_intersect_epsilon=0);
template<class T> bool Fuzzy_Intersects(RAY<VECTOR<T,2> >& ray,const ORIENTED_BOX<VECTOR<T,2> >& box,const T segment_intersect_epsilon=0);
//#####################################################################
};
};
#endif
