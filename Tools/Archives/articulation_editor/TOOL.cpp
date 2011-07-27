//#####################################################################
// Copyright 2005, Jared Go, Ranjitha Kumar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TOOL
//##################################################################### 

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MOUSE_HANDLER.h>
#include "TOOL.h"

using namespace PhysBAM;


//#####################################################################
// Constructor
//#####################################################################
template<class T> TOOL<T>::
TOOL(OPENGL_WORLD& world)
        :opengl_world(world),mousestate(NONE)
{
}
//#####################################################################
// Perform a hit test
//#####################################################################
template<class T> void TOOL<T>::
HitTest(ARRAY<OPENGL_SELECTION*> &list,int x,int y)
{
    const int buff_size=512; 
    GLuint selectBuf[buff_size];
    GLint viewport[4];
    GLint hits;

    glGetIntegerv(GL_VIEWPORT, viewport);
    glSelectBuffer(buff_size, selectBuf);
    (void) glRenderMode (GL_SELECT);
    
    glInitNames();
    glPushName(0);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();             
    glLoadIdentity();               
    gluPickMatrix((GLfloat) x, (GLfloat) (viewport[3]-y), 5.0, 5.0, viewport);
    opengl_world.Render_World(true,false);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glFlush();

    hits = glRenderMode (GL_RENDER);

    // Get selections from hit and place into list
    opengl_world.Get_Selections(list, hits, selectBuf);
}
//#####################################################################
// Handle_Click
//#####################################################################
template<class T> void TOOL<T>::
Handle_Click(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed)
{
    // Selection array
    ARRAY<OPENGL_SELECTION*> selections;

    // Depending on state machine
    switch (mousestate)
    {
    case NONE:

        // This represent no mouse buttons down state
        if (state == GLUT_UP) 
        {
            std::cout<<"Got UP state when doing nothing.. how?\n";
            return;
        }

        // Button is going down, see what we hit        
        HitTest(selections, x, y);

        // Pick the selection that we want, set that as context
        if (selections.m > 0) 
        {            
            mousestate = INTERACTING;            

            // Store what we are interacting with
            if (selections(1)->Actual_Type() == 0x1000) 
            {
                sel = selections(1);
            } 
            else if (selections(1)->Actual_Type() == 0x1001) 
            {
                sel = selections(1);
            } 
            else if (selections(1)->Actual_Type() == 0x1002) 
            {
                sel = selections(1);
            }
            opengl_world.Redisplay();
        } 
        else 
        {            
            mousestate = CHANGING_VIEW;
            opengl_world.Handle_Click_Main(button, state, x, y, ctrl_pressed, shift_pressed);
        }

        break;

    case CHANGING_VIEW:
        
        // Only respond to mouse ups
        if (state == GLUT_UP) 
        {
            opengl_world.Handle_Click_Main(button, state, x, y, ctrl_pressed, shift_pressed);
            mousestate = NONE;
        }
        break;

    case INTERACTING:

        // Move to state when released
        if (state == GLUT_UP) 
        {
            mousestate = NONE;
            sel = 0;
        }

        break;
    }

    oldmousex=x;
    oldmousey=y;
}
//#####################################################################
// Handle_Drag
//#####################################################################
template<class T> void TOOL<T>::
Handle_Drag(int x,int y)
{
    //int dx = x-oldmousex;
    //int dy = y-oldmousey;

    switch (mousestate)
    {
    case NONE:
        break;

    case CHANGING_VIEW:
        opengl_world.Handle_Drag_Main(x, y);
        break;

    case INTERACTING:

        if (sel != NULL && (sel->Actual_Type() == 0x1001 || sel->Actual_Type() == 0x1002)) 
        {                 
            OPENGL_COMPONENT_MUSCLE_3D<T> *muscle = (OPENGL_COMPONENT_MUSCLE_3D<T>*)sel->object;
            
            // Get via point in world space
            int selpt = muscle->targpoint;
            VECTOR<T,3> world = muscle->frames(selpt) * muscle->points(selpt);

            // Get opengl transform matrices
            GLint viewport[4];
            GLdouble mvmatrix[16] = {   1, 0, 0, 0,
                                        0, 1, 0, 0,
                                        0, 0, 1, 0,
                                        0, 0, 0, 1 }, projmatrix[16];
            GLint realy;             //  OpenGL y coordinate position
            GLdouble wx, wy, wz;  //  first point returned world x, y, z coords
            GLdouble win_x, win_y, win_z; 
                       
            // User Selects Guidance Points that are part of the system.            
            glGetIntegerv (GL_VIEWPORT, viewport);
            glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
            glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

            // Swap to get true y
            realy = viewport[3] - (GLint) y - 1;
            gluProject((GLdouble) world.x, (GLdouble) world.y, (GLdouble) world.z,
            mvmatrix, projmatrix, viewport, &win_x, &win_y, &win_z);

            // Project back along screen plane
            gluUnProject((GLdouble) x, (GLdouble) realy, win_z,
            mvmatrix, projmatrix, viewport, &wx, &wy, &wz);

            // Now update pos            
            std::cout<<"World is "<<world<<"\n";
            muscle->points(selpt) = muscle->frames(selpt).Inverse() * VECTOR<T,3>(wx,wy,wz);            
        }

        break;
    }

    oldmousex=x;
    oldmousey=y;
}
//#####################################################################
template class TOOL<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TOOL<double>;
#endif
