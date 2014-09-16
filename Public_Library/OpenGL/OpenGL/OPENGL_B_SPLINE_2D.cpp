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
    :curve(curve_input),bezier_version(),bezier_opengl(bezier_version,color_input),color(color_input)
, vertex_color(OPENGL_COLOR::Green(0.9)),vertex_position_color(OPENGL_COLOR::Magenta()),velocity_color(OPENGL_COLOR::Cyan()),
    draw_vertices(false),draw_velocities(false),velocity_scale(0.025)
,current_selection(0)
{
    ARRAY<TV> knot_sites(curve.knots.m-6);
    for(int i=3;i<curve.knots.m-3;i++)
        knot_sites(i-3)=curve.Evaluate(curve.knots(i));
//    LOG::printf("knot sites: %P\n",knot_sites);
    Smooth_Fit<TV>(bezier_version,knot_sites);
//    Smooth_Fit<TV>(bezier_version,curve.particles.X);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,int d> void OPENGL_B_SPLINE_2D<T,d>::
Display() const
{
    ARRAY<TV> knot_sites(curve.knots.m-6);
    for(int i=3;i<curve.knots.m-3;i++)
        knot_sites(i-3)=curve.Evaluate(curve.knots(i));
    Smooth_Fit<TV>(bezier_version,knot_sites);
    bezier_opengl.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,int d> RANGE<VECTOR<T,3> > OPENGL_B_SPLINE_2D<T,d>::
Bounding_Box() const
{
    return bezier_opengl.Bounding_Box();
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,int d> OPENGL_SELECTION<T>* OPENGL_B_SPLINE_2D<T,d>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    return bezier_opengl.Get_Selection(buffer,buffer_size);
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,int d> void OPENGL_B_SPLINE_2D<T,d>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    bezier_opengl.Highlight_Selection(selection);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,int d> void OPENGL_B_SPLINE_2D<T,d>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const
{
    bezier_opengl.Print_Selection_Info(output_stream,selection);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,int d> void OPENGL_B_SPLINE_2D<T,d>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection,MATRIX<T,3>* transform) const
{
    bezier_opengl.Print_Selection_Info(output_stream,selection,transform);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,int d> void OPENGL_B_SPLINE_2D<T,d>::
Clear_Highlight()
{
    bezier_opengl.Clear_Highlight();
}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T,int d> OPENGL_SELECTION<T>* OPENGL_B_SPLINE_2D<T,d>::
Get_Vertex_Selection(int index)
{
    return bezier_opengl.Get_Vertex_Selection(index);
}
//#####################################################################
// Function Get_Segment_Selection
//#####################################################################
template<class T,int d> OPENGL_SELECTION<T>* OPENGL_B_SPLINE_2D<T,d>::
Get_Segment_Selection(int index)
{
    return bezier_opengl.Get_Segment_Selection(index);
}
//#####################################################################
// Function Draw_Highlighted_Spline
//#####################################################################
template<class T,int d> void OPENGL_B_SPLINE_2D<T,d>::
Draw_Highlighted_Spline(int id) const
{
    bezier_opengl.Draw_Highlighted_Spline(id);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,int d> RANGE<VECTOR<T,3> > OPENGL_SELECTION_B_SPLINE_VERTEX_2D<T,d>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const BEZIER_SPLINE<VECTOR<T,2>,d> &curve=*(((OPENGL_B_SPLINE_2D<T,d> *)object)->bezier_version);
    return object->World_Space_Box(RANGE<VECTOR<T,2> >(curve.particles.X(index)));
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,int d> RANGE<VECTOR<T,3> > OPENGL_SELECTION_B_SPLINE_SEGMENT_2D<T,d>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const BEZIER_SPLINE<VECTOR<T,2>,d> &curve=*(((OPENGL_B_SPLINE_2D<T,d> *)object)->bezier_version);
    return object->World_Space_Box(RANGE<VECTOR<T,2> >::Bounding_Box(curve.particles.X.Subset(curve.control_points(index))));
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_B_SPLINE_2D<float,3>;
template class OPENGL_B_SPLINE_2D<double,3>;
}
