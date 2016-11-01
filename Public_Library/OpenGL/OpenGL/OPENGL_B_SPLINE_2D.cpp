//#####################################################################
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <OpenGL/OpenGL/OPENGL_B_SPLINE_2D.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> OPENGL_B_SPLINE_2D<T,d>::
OPENGL_B_SPLINE_2D(STREAM_TYPE stream_type,const B_SPLINE<TV,d>& curve_input,const OPENGL_COLOR &color_input)
    :BASE(stream_type,bezier_version,color_input),bezier_version(),curve(curve_input)
{
    Fill_Bezier<TV>(bezier_version,curve);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,int d> void OPENGL_B_SPLINE_2D<T,d>::
Display() const
{
    Fill_Bezier<TV>(bezier_version,curve);
    BASE::Display();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_B_SPLINE_2D<float,3>;
template class OPENGL_B_SPLINE_2D<double,3>;
}
