//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Sergey Koltakov, Neil Molino, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_SEGMENTED_CURVE_2D<T>::
OPENGL_SEGMENTED_CURVE_2D(STREAM_TYPE stream_type,const SEGMENTED_CURVE_2D<T>& curve_input,const OPENGL_COLOR &color_input)
    :OPENGL_OBJECT<T>(stream_type),curve(curve_input),color(color_input),
    vertex_color(OPENGL_COLOR::Green(0.9)),vertex_position_color(OPENGL_COLOR::Magenta()),velocity_color(OPENGL_COLOR::Cyan()),
    draw_vertices(false),draw_velocities(false),velocity_scale(0.025),selected_vertex(-1),selected_segment(-1)

{}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
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
    if(mode==GL_SELECT)
    {
        glPushName(0);
        Draw_Vertices_For_Selection();
        glLoadName(1);
        Draw_Segments_For_Selection();
        glPopName();
    }
    else
    {
        OpenGL_Begin(GL_LINES);
        for(int t=0;t<curve.mesh.elements.m;t++){
            int i=curve.mesh.elements(t)(0),j=curve.mesh.elements(t)(1);
            OpenGL_Line(curve.particles.X(i),curve.particles.X(j));}
        OpenGL_End();

        if(selected_vertex>=0)
            OPENGL_SELECTION::Draw_Highlighted_Vertex(curve.particles.X(selected_vertex),selected_vertex);
        else if(selected_segment>=0){
            int node1,node2;
            curve.mesh.elements(selected_segment).Get(node1,node2);
            OPENGL_SELECTION::Draw_Highlighted_Segment(curve.particles.X(node1),curve.particles.X(node2),selected_segment);}

        if(draw_vertices) {
            vertex_color.Send_To_GL_Pipeline();
            glPointSize(5.0f);
            OpenGL_Begin(GL_POINTS);
            for(int t=0;t<curve.particles.Size();t++)
                OpenGL_Vertex(curve.particles.X(t));
            OpenGL_End();}

        if(draw_velocities && curve.particles.store_velocity){
            velocity_color.Send_To_GL_Pipeline();
            OpenGL_Begin(GL_LINES);
            for(int t=0;t<curve.particles.Size();t++)
                OPENGL_SHAPES::Draw_Arrow(curve.particles.X(t),curve.particles.X(t)+velocity_scale*curve.particles.V(t));
            OpenGL_End();}
    }

    glPopAttrib();

    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SEGMENTED_CURVE_2D<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<TV>::Bounding_Box(curve.particles.X));
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_SEGMENTED_CURVE_2D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    PHYSBAM_ASSERT(indices.m==2);
    const static int priority[]={79,78};
    return priority[indices(0)];
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_SEGMENTED_CURVE_2D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(indices(0)==0) selected_vertex=indices(1);
    else if(indices(0)==1) selected_segment=indices(1);
    else return false;
    return true;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    Print_Selection_Info(output_stream,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Print_Selection_Info(std::ostream &output_stream,MATRIX<T,3>* transform) const
{
    if(selected_vertex>=0){
        output_stream<<"Vertex "<<selected_vertex<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(selected_vertex))<<std::endl;}
        curve.particles.Print(output_stream,selected_vertex);}
    else if(selected_segment>=0){
        int node1,node2;curve.mesh.elements(selected_segment).Get(node1,node2);
        output_stream<<"Segment "<<selected_segment<<" ("<<node1<<", "<<node2<<")"<<std::endl;
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node1<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(node1))<<std::endl;}
        curve.particles.Print(output_stream,node1);
        output_stream<<std::endl;
        output_stream<<"Vertex "<<node2<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(node2))<<std::endl;}
        curve.particles.Print(output_stream,node2);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Clear_Selection()
{
    selected_vertex=-1;
    selected_segment=-1;
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T>
void OPENGL_SEGMENTED_CURVE_2D<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION::Draw_Vertices_For_Selection(curve.mesh,curve.particles);
}
//#####################################################################
// Function Draw_Segments_For_Selection
//#####################################################################
template<class T>
void OPENGL_SEGMENTED_CURVE_2D<T>::
Draw_Segments_For_Selection() const
{
    glPushAttrib(GL_LINE_BIT);
    glLineWidth(OPENGL_PREFERENCES::selection_line_width);
    glPushName(0);
    for(int i=0;i<curve.mesh.elements.m;i++){
        int node1,node2;curve.mesh.elements(i).Get(node1,node2);
        glLoadName(i);
        OpenGL_Begin(GL_LINES);
        OpenGL_Line(curve.particles.X(node1),curve.particles.X(node2));
        OpenGL_End();
    }
    glPopName();
    glPopAttrib();
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SEGMENTED_CURVE_2D<T>::
Selection_Bounding_Box() const
{
    if(selected_vertex>=0) return World_Space_Box(RANGE<TV>(curve.particles.X(selected_vertex)));
    if(selected_segment>=0) return World_Space_Box(RANGE<TV>::Bounding_Box(curve.particles.X.Subset(curve.mesh.elements(selected_segment))));
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_SEGMENTED_CURVE_2D<float>;
template class OPENGL_SEGMENTED_CURVE_2D<double>;
}
