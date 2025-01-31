//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <OpenGL/OpenGL/OPENGL_TEXTURED_RECT.h>
#include <iostream>

using namespace PhysBAM;
using namespace std;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_TEXTURED_RECT<T>::
OPENGL_TEXTURED_RECT()
    :texture(0)
{
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_TEXTURED_RECT<T>::
Display() const
{
    if(!texture) return;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ALL_ATTRIB_BITS);

    // Force it into fill mode (e.g. it's in wireframe mode)
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glBindTexture(GL_TEXTURE_2D, texture->id);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);

    OpenGL_Begin(GL_TRIANGLE_STRIP);
    OpenGL_Texture(VECTOR<float,2>(texture->min_s,texture->min_t));
    OpenGL_Vertex(VECTOR<double,2>(-0.5*width,-0.5*height));
    OpenGL_Texture(VECTOR<float,2>(texture->min_s,texture->max_t));
    OpenGL_Vertex(VECTOR<double,2>(-0.5*width,0.5*height));
    OpenGL_Texture(VECTOR<float,2>(texture->max_s,texture->min_t));
    OpenGL_Vertex(VECTOR<double,2>(0.5*width,-0.5*height));
    OpenGL_Texture(VECTOR<float,2>(texture->max_s,texture->max_t));
    OpenGL_Vertex(VECTOR<double,2>(0.5*width,0.5*height));
    OpenGL_End();

    glDisable(GL_TEXTURE_2D);
    glEnable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);

    glPopAttrib();

    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_TEXTURED_RECT<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<T,2> >(VECTOR<T,2>(-(T)0.5*width,-(T)0.5*height),VECTOR<T,2>((T)0.5*width,(T)0.5*height)));
}
namespace PhysBAM{
template class OPENGL_TEXTURED_RECT<double>;
template class OPENGL_TEXTURED_RECT<float>;
}
