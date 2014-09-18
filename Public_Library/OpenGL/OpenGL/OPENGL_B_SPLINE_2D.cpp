//#####################################################################
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Matrices/MATRIX_3X3.h>
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
OPENGL_B_SPLINE_2D(const B_SPLINE<TV,d>& curve_input,const OPENGL_COLOR &color_input)
    :BASE(bezier_version,color_input),bezier_version(),curve(curve_input)
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
// Function Bounding_Box
//#####################################################################
template<class T,int d> RANGE<VECTOR<T,3> > OPENGL_SELECTION_B_SPLINE_VERTEX_2D<T,d>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const BEZIER_SPLINE<VECTOR<T,2>,d> &curve=((OPENGL_B_SPLINE_2D<T,d> *)object)->bezier_version;
    return object->World_Space_Box(RANGE<VECTOR<T,2> >(curve.particles.X(index)));
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,int d> RANGE<VECTOR<T,3> > OPENGL_SELECTION_B_SPLINE_SEGMENT_2D<T,d>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const BEZIER_SPLINE<VECTOR<T,2>,d> &curve=((OPENGL_B_SPLINE_2D<T,d> *)object)->bezier_version;
    return object->World_Space_Box(RANGE<VECTOR<T,2> >::Bounding_Box(curve.particles.X.Subset(curve.control_points(index))));
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_B_SPLINE_2D<float,3>;
template class OPENGL_B_SPLINE_2D<double,3>;
}
