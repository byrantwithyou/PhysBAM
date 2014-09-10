//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Sergey Koltakov, Neil Molino, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_SEGMENTED_CURVE_2D<T>::
OPENGL_SEGMENTED_CURVE_2D(const SEGMENTED_CURVE_2D<T>& curve_input,const OPENGL_COLOR &color_input)
    :curve(curve_input),color(color_input),
    vertex_color(OPENGL_COLOR::Green(0.9)),vertex_position_color(OPENGL_COLOR::Magenta()),velocity_color(OPENGL_COLOR::Cyan()),
    draw_vertices(false),draw_velocities(false),velocity_scale(0.025),current_selection(0)
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
        glPushName(1);
        Draw_Vertices_For_Selection();
        glLoadName(2);
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

        if(current_selection) {
            if(current_selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_VERTEX_2D) {
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T> *)current_selection)->index;
                OPENGL_SELECTION<T>::Draw_Highlighted_Vertex(curve.particles.X(index),index);}
            else if(current_selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_2D) {
                int index=((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T> *)current_selection)->index;
                int node1,node2;curve.mesh.elements(index).Get(node1,node2);
                OPENGL_SELECTION<T>::Draw_Highlighted_Segment(curve.particles.X(node1),curve.particles.X(node2),index);}}
    }

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
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_SEGMENTED_CURVE_2D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION<T>* selection = 0;
    if(buffer_size==2)
    {
        if(buffer[0]==1) selection = Get_Vertex_Selection(buffer[1]);
        else if(buffer[0]==2) selection = Get_Segment_Selection(buffer[1]);
    }
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    delete current_selection; current_selection = 0;
    // Make a copy of selection
    if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_VERTEX_2D)
        current_selection=new OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T>(this,((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T> *)selection)->index);
    else if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_2D)
        current_selection=new OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T>(this,((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T> *)selection)->index);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const
{
    Print_Selection_Info(output_stream,selection,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_SEGMENTED_CURVE_2D<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection,MATRIX<T,3>* transform) const
{
    if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_VERTEX_2D) {
        int index=((OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T> *)selection)->index;
        output_stream<<"Vertex "<<index<<std::endl;
        if(transform){output_stream<<"WORLD Position "<<transform->Homogeneous_Times(curve.particles.X(index))<<std::endl;}
        curve.particles.Print(output_stream,index);}
    else if(selection->type==OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_2D) {
        int index=((OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T> *)selection)->index;
        int node1,node2;curve.mesh.elements(index).Get(node1,node2);
        output_stream<<"Segment "<<index<<" ("<<node1<<", "<<node2<<")"<<std::endl;
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
Clear_Highlight()
{
    delete current_selection; current_selection = 0;
}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_SEGMENTED_CURVE_2D<T>::
Get_Vertex_Selection(int index)
{
    return new OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T>(this, index);
}
//#####################################################################
// Function Get_Segment_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_SEGMENTED_CURVE_2D<T>::
Get_Segment_Selection(int index)
{
    return new OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T>(this, index);
}
//#####################################################################
// Function Draw_Vertices_For_Selection
//#####################################################################
template<class T>
void OPENGL_SEGMENTED_CURVE_2D<T>::
Draw_Vertices_For_Selection() const
{
    OPENGL_SELECTION<T>::Draw_Vertices_For_Selection(curve.mesh,curve.particles);
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
// Function Bounding_Box
//#####################################################################
template<class T>
RANGE<VECTOR<T,3> > OPENGL_SELECTION_SEGMENTED_CURVE_VERTEX_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const SEGMENTED_CURVE_2D<T> &curve=((OPENGL_SEGMENTED_CURVE_2D<T> *)object)->curve;
    return object->World_Space_Box(RANGE<VECTOR<T,2> >(curve.particles.X(index)));
}
template<class T>
//#####################################################################
// Function Bounding_Box
//#####################################################################
RANGE<VECTOR<T,3> > OPENGL_SELECTION_SEGMENTED_CURVE_SEGMENT_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const SEGMENTED_CURVE_2D<T> &curve=((OPENGL_SEGMENTED_CURVE_2D<T> *)object)->curve;
    return object->World_Space_Box(RANGE<VECTOR<T,2> >::Bounding_Box(curve.particles.X.Subset(curve.mesh.elements(index))));
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_SEGMENTED_CURVE_2D<float>;
template class OPENGL_SEGMENTED_CURVE_2D<double>;
}
