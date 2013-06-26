//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINT_SIMPLICES_1D.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_POINT_SIMPLICES_1D<T>::
OPENGL_POINT_SIMPLICES_1D(const POINT_SIMPLICES_1D<T>& simplices_input,const OPENGL_COLOR &color_input)
    :simplices(simplices_input),color(color_input),color_gray(color_input.Grayscale()),
    vertex_color(OPENGL_COLOR::Green(0.9)),segment_color(OPENGL_COLOR::Blue(0.9)),vertex_position_color(OPENGL_COLOR::Magenta()),
    draw_vertices(false),draw_vertex_positions(false)
{}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_POINT_SIMPLICES_1D<T>::
Display(const int in_color) const
{
    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    if(in_color) color.Send_To_GL_Pipeline();
    else color_gray.Send_To_GL_Pipeline();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    if(draw_vertices){
        segment_color.Send_To_GL_Pipeline();
        glPointSize(5.0f);
        vertices.Resize(0);
        for(int t=0;t<simplices.particles.Size();t++){
            OpenGL_Vertex(simplices.particles.X(t),vertices);}
        OpenGL_Draw_Arrays(GL_LINES,2,vertices);
        vertex_color.Send_To_GL_Pipeline();
        glPointSize(10.0f);
        vertices.Resize(0);
        for(int t=0;t<simplices.particles.Size();t++){
            OpenGL_Vertex(simplices.particles.X(t),vertices);}
        OpenGL_Draw_Arrays(GL_POINTS,2,vertices);}

#ifndef USE_OPENGLES
    if (draw_vertex_positions) {
        vertex_position_color.Send_To_GL_Pipeline();
        for(int t=0;t<simplices.particles.Size();t++){
            VECTOR<float,3> world_space_pos=World_Space_Point(VECTOR<float,2>(simplices.particles.X(t)));
            OpenGL_String(simplices.particles.X(t),STRING_UTILITIES::string_sprintf("<%f>",world_space_pos.x));}}
#endif

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_POINT_SIMPLICES_1D<T>::
Bounding_Box() const
{
    float xmin,xmax;
    xmin=xmax=simplices.particles.X(1).x;
    for(int i=0;i<simplices.particles.Size();i++){
        xmin=min(xmin,(float)simplices.particles.X(i).x);xmax=max(xmax,(float)simplices.particles.X(i).x);}
    return World_Space_Box(RANGE<VECTOR<float,3> >(VECTOR<float,3>(xmin,0,0),VECTOR<float,3>(xmax,0,0)));
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_POINT_SIMPLICES_1D<float>;
template class OPENGL_POINT_SIMPLICES_1D<double>;
}
