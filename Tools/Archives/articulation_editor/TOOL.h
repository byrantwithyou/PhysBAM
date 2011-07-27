//#####################################################################
// Copyright 2005, Jared Go, Ranjitha Kumar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TOOL
//##################################################################### 
#ifndef __TOOL__
#define __TOOL__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MOUSE_HANDLER.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include "OPENGL_COMPONENT_MUSCLE_3D.h"

namespace PhysBAM
{

template<class T>
class TOOL:public OPENGL_MOUSE_HANDLER
{
public:
    OPENGL_WORLD& opengl_world;
    
    // Our "selection" - muscles for now, we want an interface moveable
    // does this live 
    OPENGL_SELECTION* sel;    
 
    // Mouse FSM
    enum STATE { NONE, CHANGING_VIEW, INTERACTING };
    STATE mousestate;

    // Mouse
    int oldmousex,oldmousey;    
    
    //Members
    TOOL(OPENGL_WORLD&);
    virtual void Handle_Click(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed);
    virtual void Handle_Drag(int x,int y);

private:
    
    void HitTest(ARRAY<OPENGL_SELECTION*> &list, int x, int y);

};

}
#endif
