//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL/OPENGL_WINDOW_GLUT.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>

using namespace PhysBAM;

//#####################################################################
// OPENGL_WINDOW_GLUT
//#####################################################################
template<class T> OPENGL_WINDOW_GLUT<T>::
OPENGL_WINDOW_GLUT(OPENGL_WORLD<T>& opengl_world_input,const std::string& window_title_input,const int width_input,const int height_input)
    :OPENGL_WINDOW<T>(opengl_world_input),width(width_input),height(height_input)
{
    if(Single()) PHYSBAM_FATAL_ERROR("Only one glut context allowed");
    Single()=this;

    static int argc=1;static const char *(argv[1]);argv[0]="Visualization";
    glutInit(&argc,(char**)argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_ALPHA);
    glutInitWindowSize(width,height);
    main_window=glutCreateWindow(window_title_input.c_str());
    glutDisplayFunc(Handle_Display_Glut);
    glutReshapeFunc(Handle_Reshape_Glut);
    glutKeyboardFunc(Handle_Keypress_Glut);
    glutSpecialFunc(Handle_Special_Keypress_Glut);
    glutMouseFunc(Handle_Click_Glut);
    glutMotionFunc(Handle_Drag_Glut);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_WINDOW_GLUT<T>::
~OPENGL_WINDOW_GLUT()
{Single()=0;}
//#####################################################################
// Setup_Idle
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Setup_Idle(const bool use)
{
    glutIdleFunc(use?Handle_Idle_Glut:0);
}
//#####################################################################
// Setup_Timer
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Setup_Timer(const float wait_milliseconds)
{
    glutTimerFunc((int)(wait_milliseconds*1000)+1,Handle_Timer_Glut,0);
}
//#####################################################################
// Redisplay
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Redisplay()
{
    glutPostRedisplay();
}
//#####################################################################
// Main_Loop
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Main_Loop()
{
    glutMainLoop();
}
//#####################################################################
// Request_Resize
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Request_Resize(const int width,const int height)
{
    glutReshapeWindow(width,height);
}
//#####################################################################
// Request_Move
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Request_Move(const int x,const int y)
{
    glutPositionWindow(x,y);
}
//#####################################################################
// Width
//#####################################################################
template<class T> int OPENGL_WINDOW_GLUT<T>::
Width() const
{
    return width;
}
//#####################################################################
// Height
//#####################################################################
template<class T> int OPENGL_WINDOW_GLUT<T>::
Height() const
{
    return height;
}
//#####################################################################
// Handle_Display_Glut
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Handle_Display_Glut()
{
    Single()->opengl_world.Render_World(false); // render, no selection
    glutSwapBuffers();
}
//#####################################################################
// Handle_Reshape_Glut
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Handle_Reshape_Glut(int w,int h)
{
    Single()->width=w;Single()->height=h;
    Single()->opengl_world.Handle_Reshape_Main();
}

//#####################################################################
// Handle_Special_Keypress_Glut
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Handle_Special_Keypress_Glut(int key,int x,int y)
{
    Single()->opengl_world.Handle_Keypress_Main(OPENGL_KEY::From_Glut_Special_Key(key,(glutGetModifiers()&GLUT_ACTIVE_CTRL)!=0,(glutGetModifiers()&GLUT_ACTIVE_ALT)!=0),x,y);
}
//#####################################################################
// Handle_Keypress_Glut
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Handle_Keypress_Glut(unsigned char key,int x,int y)
{
    if(Single()->opengl_world.prompt_mode) Single()->opengl_world.Handle_Keypress_Prompt(key);
    else Single()->opengl_world.Handle_Keypress_Main(OPENGL_KEY::From_Glut_Key(key,(glutGetModifiers()&GLUT_ACTIVE_CTRL)!=0,(glutGetModifiers()&GLUT_ACTIVE_ALT)!=0),x,y);
}
//#####################################################################
// Handle_Click_Glut
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Handle_Click_Glut(int button,int state,int x,int y)
{
    bool ctrl_pressed=(glutGetModifiers() & GLUT_ACTIVE_CTRL)!=0;
    bool shift_pressed=glutGetModifiers() & GLUT_ACTIVE_SHIFT;
    Single()->opengl_world.Handle_Click_Main(button,state,x,y,ctrl_pressed,shift_pressed);
}
//#####################################################################
// Handle_Drag_Glut
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Handle_Drag_Glut(int x,int y)
{
    Single()->opengl_world.Handle_Drag_Main(x,y);
}
//#####################################################################
// Handle_Idle_Glut
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Handle_Idle_Glut()
{
    Single()->opengl_world.Handle_Idle();
}
//#####################################################################
// Handle_Timer_Glut
//#####################################################################
template<class T> void OPENGL_WINDOW_GLUT<T>::
Handle_Timer_Glut(int value)
{
    Single()->opengl_world.Handle_Timer();
}
//#####################################################################
// Function Single
//#####################################################################
template<class T> OPENGL_WINDOW_GLUT<T>*& OPENGL_WINDOW_GLUT<T>::
Single()
{
    static OPENGL_WINDOW_GLUT<T>* single=0;
    return single;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_WINDOW_GLUT<float>;
template class OPENGL_WINDOW_GLUT<double>;
}
