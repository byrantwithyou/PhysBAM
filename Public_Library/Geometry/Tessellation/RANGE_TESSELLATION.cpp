//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <Core/Math_Tools/RANGE.h>
#include <Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(
    const RANGE<VECTOR<T,3> >& box,const VECTOR<int,3>& res)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    surface->Initialize_Box_Mesh_And_Particles(res,box);
    return surface;
}
//#####################################################################
}
template TRIANGULATED_SURFACE<float>* TESSELLATION::Generate_Triangles(
    const RANGE<VECTOR<float,3> >&,const VECTOR<int,3>&);
template TRIANGULATED_SURFACE<double>* TESSELLATION::Generate_Triangles(
    const RANGE<VECTOR<double,3> >&,const VECTOR<int,3>&);
}
