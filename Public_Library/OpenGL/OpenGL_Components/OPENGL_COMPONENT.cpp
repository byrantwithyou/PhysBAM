//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
OPENGL_COMPONENT::
OPENGL_COMPONENT(const std::string &name)
 : frame(0), draw(true), is_animation(false), component_name(name)
{
}
//#####################################################################
// Destructor
//#####################################################################
OPENGL_COMPONENT::
~OPENGL_COMPONENT()
{
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
bool OPENGL_COMPONENT::
Valid_Frame(int frame_input) const
{
    return false;
}
//#####################################################################
// Function Is_Up_To_Date
//#####################################################################
bool OPENGL_COMPONENT::
Is_Up_To_Date(int frame) const
{
    return true;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
void OPENGL_COMPONENT::
Set_Frame(int frame_input)
{
    frame = frame_input;
}
//#####################################################################
// Function Set_Draw
//#####################################################################
void OPENGL_COMPONENT::
Set_Draw(bool draw_input)
{
    draw = draw_input;
}
//#####################################################################
// Function Draw_All_Objects
//#####################################################################
void OPENGL_COMPONENT::
Draw_All_Objects()
{
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
bool OPENGL_COMPONENT::
Use_Bounding_Box() const
{
    return draw;
}
//#####################################################################
// Function Next_Frame
//#####################################################################
void OPENGL_COMPONENT::
Next_Frame()
{
    Set_Frame(frame+1);
}
//#####################################################################
// Function Prev_Frame
//#####################################################################
void OPENGL_COMPONENT::
Prev_Frame()
{
    Set_Frame(frame-1);
}
