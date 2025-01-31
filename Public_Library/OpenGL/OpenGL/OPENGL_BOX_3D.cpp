//#####################################################################
// Copyright 2003, 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <OpenGL/OpenGL/OPENGL_BOX_3D.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>

using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_BOX_3D<T>::
Display() const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    color.Send_To_GL_Pipeline();
    glLineWidth(3.0f);

    VECTOR<T,3> x_vector(box.max_corner.x-box.min_corner.x,0,0), y_vector(0,box.max_corner.y-box.min_corner.y,0), z_vector(0,0,box.max_corner.z-box.min_corner.z);
    VECTOR<T,3> left_corner(box.min_corner.x,box.min_corner.y,box.min_corner.z);
    VECTOR<T,3> right_corner=left_corner+x_vector;
    OpenGL_Begin(GL_LINE_LOOP);
    OpenGL_Vertex(left_corner);
    OpenGL_Vertex(left_corner+y_vector); 
    OpenGL_Vertex(left_corner+y_vector+z_vector);
    OpenGL_Vertex(left_corner+z_vector);
    OpenGL_End();
    OpenGL_Begin(GL_LINE_LOOP);
    OpenGL_Vertex(right_corner);
    OpenGL_Vertex(right_corner+y_vector); 
    OpenGL_Vertex(right_corner+y_vector+z_vector);
    OpenGL_Vertex(right_corner+z_vector);
    OpenGL_End();
    OpenGL_Begin(GL_LINES);
    OpenGL_Line(left_corner,right_corner);
    OpenGL_Line(left_corner+y_vector,right_corner+y_vector);
    OpenGL_Line(left_corner+y_vector+z_vector,right_corner+y_vector+z_vector);
    OpenGL_Line(left_corner+z_vector,right_corner+z_vector);
    OpenGL_End();

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_BOX_3D<T>::
Bounding_Box() const
{
    return World_Space_Box(box);
}
namespace PhysBAM{
template class OPENGL_BOX_3D<float>;
template class OPENGL_BOX_3D<double>;
}
