//#####################################################################
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <OpenGL/OpenGL/OPENGL_BEZIER_SPLINE_2D.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> OPENGL_BEZIER_SPLINE_2D<T,d>::
OPENGL_BEZIER_SPLINE_2D(STREAM_TYPE stream_type,const BEZIER_SPLINE<TV,d>& curve_input,const OPENGL_COLOR &color_input)
    :OPENGL_OBJECT<T>(stream_type),curve(curve_input),color(color_input),
    vertex_color(OPENGL_COLOR::Green(0.9)),vertex_position_color(OPENGL_COLOR::Magenta()),velocity_color(OPENGL_COLOR::Cyan()),
    draw_vertices(false),draw_velocities(false),velocity_scale(0.025),selected_vertex(-1),selected_segment(-1)
{}
//#####################################################################
// Function Display
//#####################################################################
template<class T,int d> void OPENGL_BEZIER_SPLINE_2D<T,d>::
Display() const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    color.Send_To_GL_Pipeline();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    GLint mode=0;
    glGetIntegerv(GL_RENDER_MODE, &mode);
    if(mode==GL_SELECT){
        glPushName(0);
        Draw_Vertices_For_Selection();
        glLoadName(1);
        Draw_Curves_For_Selection();
        glPopName();}
    else
    {
        for(int t=0;t<curve.control_points.m;t++){
            VECTOR<TV,d+1> control_points(curve.particles.X.Subset(curve.control_points(t)));
            OpenGL_Draw_Spline(GL_LINE,20,control_points);}

        if(selected_vertex>=0) OPENGL_SELECTION::Draw_Highlighted_Vertex(curve.particles.X(selected_vertex),selected_vertex);
        if(selected_segment>=0) Draw_Highlighted_Spline(selected_segment);
    }

    if(draw_vertices){
        vertex_color.Send_To_GL_Pipeline();
        glPointSize(5.0f);
        OpenGL_Begin(GL_POINTS);
        for(int t=0;t<curve.particles.Size();t++){
            OpenGL_Vertex(curve.particles.X(t));}
        OpenGL_End();}

    if(draw_velocities && curve.particles.store_velocity){
        velocity_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int t=0;t<curve.particles.Size();t++)
            OPENGL_SHAPES::Draw_Arrow(curve.particles.X(t),curve.particles.X(t)+velocity_scale*curve.particles.V(t));
        OpenGL_End();}

    glPopAttrib();

    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,int d> RANGE<VECTOR<T,3> > OPENGL_BEZIER_SPLINE_2D<T,d>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<TV>::Bounding_Box(curve.particles.X.Subset(curve.control_points.Flattened())));
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T,int d> int OPENGL_BEZIER_SPLINE_2D<T,d>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    const static int priority[]={74,73};
    return priority[indices(0)];
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,int d> bool OPENGL_BEZIER_SPLINE_2D<T,d>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    PHYSBAM_ASSERT(indices.m==2);
    if(indices(0)==0) selected_vertex=indices(1);
    else if(indices(0)==1) selected_segment=indices(1);
    else return false;
    return true;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,int d> void OPENGL_BEZIER_SPLINE_2D<T,d>::
Print_Selection_Info(std::ostream &output_stream) const
{
    Print_Selection_Info(output_stream,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,int d> void OPENGL_BEZIER_SPLINE_2D<T,d>::
Print_Selection_Info(std::ostream &output_stream,MATRIX<T,3>* transform) const
{
    if(selected_vertex>=0){
        output_stream<<"Vertex "<<selected_vertex<<std::endl;
        if(transform)
            output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(selected_vertex))<<std::endl;
        curve.particles.Print(output_stream,selected_vertex);}
    else if(selected_segment>=0){
        CV c=curve.control_points(selected_segment);
        output_stream<<"Segment "<<selected_segment<< c <<std::endl;
        for(int i=0;i<d+1;i++){
            int node=c(i);
            output_stream<<std::endl;
            output_stream<<"Vertex "<<node<<std::endl;
            if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(node))<<std::endl;}
            curve.particles.Print(output_stream,node);}}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,int d> void OPENGL_BEZIER_SPLINE_2D<T,d>::
Clear_Selection()
{
    selected_vertex=-1;
    selected_segment=-1;
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T,int d>
void OPENGL_BEZIER_SPLINE_2D<T,d>::
Draw_Vertices_For_Selection() const
{
    glPushAttrib(GL_POINT_BIT);
    glPointSize(OPENGL_PREFERENCES::selection_point_size);
    glPushName(0);
    ARRAY<int> particles_in_mesh;
    Get_Unique(particles_in_mesh,curve.control_points.Flattened());
    ARRAY<typename OPENGL_POLICY<T>::T_GL >vertices;
    for(int i=0;i<particles_in_mesh.m;i++){const int p=particles_in_mesh(i);
        glLoadName(p);
        OpenGL_Begin(GL_POINTS);
        OpenGL_Vertex(curve.particles.X(p));
        OpenGL_End();
    }
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Curves_For_Selection
//#####################################################################
template<class T,int d> void OPENGL_BEZIER_SPLINE_2D<T,d>::
Draw_Curves_For_Selection() const
{
    glPushAttrib(GL_LINE_BIT);
    glLineWidth(OPENGL_PREFERENCES::selection_line_width);
    glPushName(0);
    ARRAY<typename OPENGL_POLICY<T>::T_GL> control_points;
    for(int t=0;t<curve.control_points.m;t++){
        glLoadName(t);
        VECTOR<TV,d+1> control_points(curve.particles.X.Subset(curve.control_points(t)));
        OpenGL_Draw_Spline(GL_LINE,20,control_points);
    }
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Draw_Highlighted_Spline
//#####################################################################
template<class T,int d> void OPENGL_BEZIER_SPLINE_2D<T,d>::
Draw_Highlighted_Spline(int id) const
{
    glPushAttrib(GL_LINE_BIT | GL_ENABLE_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);
    OPENGL_PREFERENCES::selection_highlight_color.Send_To_GL_Pipeline();
    VECTOR<TV,d+1> control_points(curve.particles.X.Subset(curve.control_points(id)));
    OpenGL_Draw_Spline(GL_LINE,20,control_points);
    if(id>=0) OpenGL_String(curve.Evaluate(id,.5),LOG::sprintf("   %d",id));
    glPopAttrib();
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T,int d> RANGE<VECTOR<T,3> > OPENGL_BEZIER_SPLINE_2D<T,d>::
Selection_Bounding_Box() const
{
    if(selected_vertex) return World_Space_Box(RANGE<TV>(curve.particles.X(selected_vertex)));
    if(selected_segment) return World_Space_Box(RANGE<TV>::Bounding_Box(curve.particles.X.Subset(curve.control_points(selected_segment))));
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_BEZIER_SPLINE_2D<float,3>;
template class OPENGL_BEZIER_SPLINE_2D<double,3>;
}
