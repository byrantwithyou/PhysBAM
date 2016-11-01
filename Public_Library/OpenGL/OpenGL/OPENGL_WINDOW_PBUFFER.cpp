//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL/OPENGL_WINDOW_PBUFFER.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#ifdef _WIN32
#include <OpenGL/OpenGL/OPENGL_WGL_PBUFFER.h>
//#elif defined(__APPLE__)
//#include <OpenGL/OpenGL/OPENGL_AGL_PBUFFER.h>
#else
#include <OpenGL/OpenGL/OPENGL_GLX_PBUFFER.h>
#endif

using namespace PhysBAM;

//#####################################################################
// OPENGL_WINDOW_PBUFFER
//#####################################################################
template<class T> OPENGL_WINDOW_PBUFFER<T>::
OPENGL_WINDOW_PBUFFER(OPENGL_WORLD<T>& opengl_world_input,const std::string& window_title_input,const int width_input,const int height_input)
    :OPENGL_WINDOW<T>(opengl_world_input),width(width_input),height(height_input)
{
    pbuffer=new OPENGL_PBUFFER;
    if(!pbuffer->Create(width,height)) PHYSBAM_FATAL_ERROR("Could not set up pbuffer");
    // need glut init to make sure we can run glut functions
    static int argc=1;static const char *(argv[1]);argv[0]="Visualization";
    glutInit(&argc,(char**)argv);
}
//#####################################################################
// ~OPENGL_WINDOW_PBUFFER
//#####################################################################
template<class T> OPENGL_WINDOW_PBUFFER<T>::
~OPENGL_WINDOW_PBUFFER()
{
    delete pbuffer;
}
//#####################################################################
// Function Handle_Idle
//#####################################################################
template<class T> void OPENGL_WINDOW_PBUFFER<T>::
Setup_Idle(const bool use)
{}
//#####################################################################
// Function Handle_Idle
//#####################################################################
template<class T> void OPENGL_WINDOW_PBUFFER<T>::
Setup_Timer(const float wait_milliseconds)
{}
//#####################################################################
// Function Handle_Idle
//#####################################################################
template<class T> void OPENGL_WINDOW_PBUFFER<T>::
Redisplay()
{}
//#####################################################################
// Function Handle_Idle
//#####################################################################
template<class T> void OPENGL_WINDOW_PBUFFER<T>::
Main_Loop()
{}
//#####################################################################
// Request_Resize
//#####################################################################
template<class T> void OPENGL_WINDOW_PBUFFER<T>::
Request_Resize(const int width,const int height)
{}
//#####################################################################
// Request_Move
//#####################################################################
template<class T> void OPENGL_WINDOW_PBUFFER<T>::
Request_Move(const int x,const int y)
{}
//#####################################################################
// Width
//#####################################################################
template<class T> int OPENGL_WINDOW_PBUFFER<T>::
Width() const
{ 
    return width;
}
//#####################################################################
// Height
//#####################################################################
template<class T> int OPENGL_WINDOW_PBUFFER<T>::
Height() const
{
    return height;
}
//#####################################################################
namespace PhysBAM
{
template class OPENGL_WINDOW_PBUFFER<double>;
template class OPENGL_WINDOW_PBUFFER<float>;
}
