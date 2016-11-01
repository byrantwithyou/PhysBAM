//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <OpenGL/OpenGL/OPENGL_POINT_SIMPLICES_1D.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_POINT_SIMPLICES_1D<T>::
OPENGL_POINT_SIMPLICES_1D(STREAM_TYPE stream_type,const POINT_SIMPLICES_1D<T>& simplices_input,const OPENGL_COLOR &color_input)
    :OPENGL_OBJECT<T>(stream_type),simplices(simplices_input),color(color_input),
    vertex_color(OPENGL_COLOR::Green(0.9)),segment_color(OPENGL_COLOR::Blue(0.9)),vertex_position_color(OPENGL_COLOR::Magenta()),
    draw_vertices(false)
{}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_POINT_SIMPLICES_1D<T>::
Display() const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    color.Send_To_GL_Pipeline();

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    if(draw_vertices){
        segment_color.Send_To_GL_Pipeline();
        glPointSize(5.0f);
        OpenGL_Begin(GL_LINES);
        for(int t=0;t<simplices.particles.Size();t++){
            OpenGL_Vertex(simplices.particles.X(t));}
        OpenGL_End();
        vertex_color.Send_To_GL_Pipeline();
        glPointSize(10.0f);
        OpenGL_Begin(GL_POINTS);
        for(int t=0;t<simplices.particles.Size();t++){
            OpenGL_Vertex(simplices.particles.X(t));}
        OpenGL_End();}

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_POINT_SIMPLICES_1D<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<T,3> >(RANGE<TV>::Bounding_Box(simplices.particles.X)));
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_POINT_SIMPLICES_1D<float>;
template class OPENGL_POINT_SIMPLICES_1D<double>;
}
