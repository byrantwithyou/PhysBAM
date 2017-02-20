//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __BOX_POLYGON_INTERSECTION__
#define __BOX_POLYGON_INTERSECTION__
#include <Core/Arrays/ARRAYS_FORWARD.h>

namespace PhysBAM{
template<class TV> class RANGE;
namespace INTERSECTION{
//#####################################################################
template<class TV> void Box_Polygon_Intersection(const RANGE<TV>& box,ARRAY<TV>& vertices,bool repeat_first_vertex);
//#####################################################################
}
}
#endif
