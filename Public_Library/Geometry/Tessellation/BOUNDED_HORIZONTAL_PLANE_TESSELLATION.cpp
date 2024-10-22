//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <Geometry/Tessellation/BOUNDED_HORIZONTAL_PLANE_TESSELLATION.h>
#include <Geometry/Tessellation/PLANE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const BOUNDED_HORIZONTAL_PLANE<VECTOR<T,3> >& plane){
    return Generate_Triangles(PLANE<T>(VECTOR<T,3>(0,1,0),VECTOR<T,3>()));
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const BOUNDED_HORIZONTAL_PLANE<VECTOR<float,3> >&);
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const BOUNDED_HORIZONTAL_PLANE<VECTOR<double,3> >&);
}
}
