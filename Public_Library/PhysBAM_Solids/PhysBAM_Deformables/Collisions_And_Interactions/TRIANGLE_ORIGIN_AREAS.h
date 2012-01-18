//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_ORIGIN_AREAS
//##################################################################### 
#ifndef __TRIANGLE_ORIGIN_AREAS__
#define __TRIANGLE_ORIGIN_AREAS__

#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SEGMENT_ORIGIN_AREAS.h>
namespace PhysBAM{
namespace ORIGIN_AREAS
{
template<class T,class TV> void Volume_From_Simplices(VOL_DATA<T,3,6>& data,const TV& X0,const TV A[6]);
}
}
#endif

