//#####################################################################
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Matrices/MATRIX_3X3.h>
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
    draw_vertices(false),draw_velocities(false),velocity_scale(0.025),current_selection(0)
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
        glPushName(1);
        Draw_Vertices_For_Selection();
        glLoadName(2);
        Draw_Curves_For_Selection();
        glPopName();}
    else
    {
        for(int t=0;t<curve.control_points.m;t++){
            VECTOR<TV,d+1> control_points(curve.particles.X.Subset(curve.control_points(t)));
            OpenGL_Draw_Spline(GL_LINE,20,control_points);}

        if(current_selection){
            if(current_selection->type==OPENGL_SELECTION<T>::BEZIER_SPLINE_VERTEX_2D){
                int index=((OPENGL_SELECTION_BEZIER_SPLINE_VERTEX_2D<T,d> *)current_selection)->index;
                OPENGL_SELECTION<T>::Draw_Highlighted_Vertex(curve.particles.X(index),index);}
            else if(current_selection->type==OPENGL_SELECTION<T>::BEZIER_SPLINE_SEGMENT_2D){
                int index=((OPENGL_SELECTION_BEZIER_SPLINE_SEGMENT_2D<T,d> *)current_selection)->index;
                Draw_Highlighted_Spline(index);}}
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
// Function Get_Selection
//#####################################################################
template<class T,int d> OPENGL_SELECTION<T>* OPENGL_BEZIER_SPLINE_2D<T,d>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION<T>* selection=0;
    if(buffer_size==2){
        if(buffer[0]==1) selection=Get_Vertex_Selection(buffer[1]);
        else if(buffer[0]==2) selection=Get_Segment_Selection(buffer[1]);}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,int d> void OPENGL_BEZIER_SPLINE_2D<T,d>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    delete current_selection; current_selection=0;
    // Make a copy of selection
    if(selection->type==OPENGL_SELECTION<T>::BEZIER_SPLINE_VERTEX_2D)
        current_selection=new OPENGL_SELECTION_BEZIER_SPLINE_VERTEX_2D<T,d>(this,((OPENGL_SELECTION_BEZIER_SPLINE_VERTEX_2D<T,d> *)selection)->index);
    else if(selection->type==OPENGL_SELECTION<T>::BEZIER_SPLINE_SEGMENT_2D)
        current_selection=new OPENGL_SELECTION_BEZIER_SPLINE_SEGMENT_2D<T,d>(this,((OPENGL_SELECTION_BEZIER_SPLINE_SEGMENT_2D<T,d> *)selection)->index);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,int d> void OPENGL_BEZIER_SPLINE_2D<T,d>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const
{
    Print_Selection_Info(output_stream,selection,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,int d> void OPENGL_BEZIER_SPLINE_2D<T,d>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection,MATRIX<T,3>* transform) const
{
    if(selection->type==OPENGL_SELECTION<T>::BEZIER_SPLINE_VERTEX_2D) {
        int index=((OPENGL_SELECTION_BEZIER_SPLINE_VERTEX_2D<T,d> *)selection)->index;
        output_stream<<"Vertex "<<index<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(index))<<std::endl;}
        curve.particles.Print(output_stream,index);}
    else if(selection->type==OPENGL_SELECTION<T>::BEZIER_SPLINE_SEGMENT_2D) {
        int index=((OPENGL_SELECTION_BEZIER_SPLINE_SEGMENT_2D<T,d> *)selection)->index;
        CV c=curve.control_points(index);
        output_stream<<"Segment "<<index<< c <<std::endl;
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
Clear_Highlight()
{
    delete current_selection;
    current_selection=0;
}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T,int d> OPENGL_SELECTION<T>* OPENGL_BEZIER_SPLINE_2D<T,d>::
Get_Vertex_Selection(int index)
{
    return new OPENGL_SELECTION_BEZIER_SPLINE_VERTEX_2D<T,d>(this, index);
}
//#####################################################################
// Function Get_Segment_Selection
//#####################################################################
template<class T,int d> OPENGL_SELECTION<T>* OPENGL_BEZIER_SPLINE_2D<T,d>::
Get_Segment_Selection(int index)
{
    return new OPENGL_SELECTION_BEZIER_SPLINE_SEGMENT_2D<T,d>(this, index);
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
    curve.control_points.Flattened().Get_Unique(particles_in_mesh);
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
    if(id>=0) OpenGL_String(curve.Evaluate(id,.5),STRING_UTILITIES::string_sprintf("   %d",id));
    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,int d> RANGE<VECTOR<T,3> > OPENGL_SELECTION_BEZIER_SPLINE_VERTEX_2D<T,d>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const BEZIER_SPLINE<VECTOR<T,2>,d> &curve=((OPENGL_BEZIER_SPLINE_2D<T,d> *)object)->curve;
    return object->World_Space_Box(RANGE<VECTOR<T,2> >(curve.particles.X(index)));
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,int d> RANGE<VECTOR<T,3> > OPENGL_SELECTION_BEZIER_SPLINE_SEGMENT_2D<T,d>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const BEZIER_SPLINE<VECTOR<T,2>,d> &curve=((OPENGL_BEZIER_SPLINE_2D<T,d> *)object)->curve;
    return object->World_Space_Box(RANGE<VECTOR<T,2> >::Bounding_Box(curve.particles.X.Subset(curve.control_points(index))));
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_BEZIER_SPLINE_2D<float,3>;
template class OPENGL_BEZIER_SPLINE_2D<double,3>;
}
