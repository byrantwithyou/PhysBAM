//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __PLANE_TESSELLATION__
#define __PLANE_TESSELLATION__
 
namespace PhysBAM{
template<class T> class PLANE;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class SEGMENTED_CURVE_2D;
template<class T> class LINE_2D;

namespace TESSELLATION{
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const PLANE<T>& plane,T tesselated_size=1000,
    int elements_per_side=1);
template<class T> TRIANGULATED_SURFACE<T>* Tessellate_Boundary(const PLANE<T>& plane,T tesselated_size=1000,
    int elements_per_side=1)
{return Generate_Triangles(plane,tesselated_size,elements_per_side);}
template<class T> SEGMENTED_CURVE_2D<T>* Tessellate_Boundary(const LINE_2D<T>& line,T tesselated_size=1000,
    int elements_per_side=1);
//#####################################################################
}
}
#endif
