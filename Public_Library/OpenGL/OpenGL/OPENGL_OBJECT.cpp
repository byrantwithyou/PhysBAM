//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Geoffrey Irving, Michael Lentine, Neil Molino, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
OPENGL_OBJECT::
OPENGL_OBJECT()
    :frame(&default_frame),selectable(false),visible(true),show_name(true),slice(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
OPENGL_OBJECT::
~OPENGL_OBJECT()
{
}
//#####################################################################
// Function Display
//#####################################################################
void OPENGL_OBJECT::
Display() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
bool OPENGL_OBJECT::
Use_Bounding_Box() const
{
    return true;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
RANGE<VECTOR<float,3> > OPENGL_OBJECT::
Bounding_Box() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Is_Transparent
//#####################################################################
bool OPENGL_OBJECT::
Is_Transparent() const
{
    return false;  // this is usually what we want
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
void OPENGL_OBJECT::
Turn_Smooth_Shading_On()
{
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
void OPENGL_OBJECT::
Turn_Smooth_Shading_Off()
{
}
//#####################################################################
// Function Get_Selection
//#####################################################################
OPENGL_SELECTION* OPENGL_OBJECT::
Get_Selection(GLuint *buffer,int buffer_size)
{
    return 0;
}
//#####################################################################
// Function Set_Selection
//#####################################################################
void OPENGL_OBJECT::
Set_Selection(OPENGL_SELECTION* selection)
{
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
void OPENGL_OBJECT::
Highlight_Selection(OPENGL_SELECTION* selection)
{
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
void OPENGL_OBJECT::
Clear_Highlight()
{
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
void OPENGL_OBJECT::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION* selection) const
{
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
RANGE<VECTOR<float,3> > OPENGL_OBJECT::
Selection_Bounding_Box(OPENGL_SELECTION* selection) const
{
    return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
OPENGL_SELECTION* OPENGL_OBJECT::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    return 0;
}
//#####################################################################
// Function Set_Slice
//#####################################################################
void OPENGL_OBJECT::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;Slice_Has_Changed();
}
//#####################################################################
// Function Slice_Has_Changed
//#####################################################################
void OPENGL_OBJECT::
Slice_Has_Changed()
{
}
