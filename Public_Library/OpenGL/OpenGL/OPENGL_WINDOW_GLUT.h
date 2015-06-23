//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OPENGL_WINDOW_GLUT__
#define __OPENGL_WINDOW_GLUT__

#include <OpenGL/OpenGL/OPENGL_WINDOW.h>
namespace PhysBAM{

template<class T>
class OPENGL_WINDOW_GLUT:public OPENGL_WINDOW<T>
{
    using OPENGL_WINDOW<T>::opengl_world;
    int main_window;
    int width,height;

//#####################################################################
public:
    OPENGL_WINDOW_GLUT(OPENGL_WORLD<T>& world_input,const std::string& window_title_input,const int width_input,const int height_input);
    virtual ~OPENGL_WINDOW_GLUT();
    void Setup_Idle(const bool use) override;
    void Setup_Timer(const float wait_milliseconds) override;
    void Redisplay() override;
    void Main_Loop() override;
    void Request_Resize(const int width,const int height) override;
    void Request_Move(const int x,const int y) override;
    int Width() const override;
    int Height() const override;
private:
    static OPENGL_WINDOW_GLUT<T>*& Single();
    static void Handle_Idle_Glut();
    static void Handle_Timer_Glut(int value);
    static void Handle_Display_Glut();
    static void Handle_Reshape_Glut(int w,int h);
    static void Handle_Special_Keypress_Glut(int key,int x,int y);
    static void Handle_Keypress_Glut(unsigned char key,int x,int y);
    static void Handle_Click_Glut(int button,int state,int x,int y);
    static void Handle_Drag_Glut(int x,int y);
//#####################################################################
};
}
#endif
