//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/FRAME.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Tessellation/ORIENTED_BOX_TESSELLATION.h>
#include <Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const ORIENTED_BOX<VECTOR<T,3> >& box)
{
    typedef VECTOR<T,3> TV;
    MATRIX<T,3> directions=box.edges;TV lengths;
    for(int i=0;i<3;i++) lengths[i]=directions.Column(i).Normalize();
    RANGE<TV> aligned_box(TV(),lengths);
    TRIANGULATED_SURFACE<T>* surface=Generate_Triangles(aligned_box);
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    FRAME<TV> frame(box.corner,ROTATION<TV>(directions));
    for(int i=0;i<particles.X.m;i++) particles.X(i)=frame*particles.X(i);
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const ORIENTED_BOX<VECTOR<float,3> >&);
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const ORIENTED_BOX<VECTOR<double,3> >&);
}
}
