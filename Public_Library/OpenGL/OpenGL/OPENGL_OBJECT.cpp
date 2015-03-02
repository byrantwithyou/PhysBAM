//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Geoffrey Irving, Michael Lentine, Neil Molino, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_OBJECT<T>::
OPENGL_OBJECT(STREAM_TYPE stream_type)
    :frame(&default_frame),selectable(false),visible(true),show_name(true),slice(0),stream_type(stream_type)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_OBJECT<T>::
~OPENGL_OBJECT()
{
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_OBJECT<T>::
Display() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_OBJECT<T>::
Use_Bounding_Box() const
{
    return true;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_OBJECT<T>::
Bounding_Box() const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Is_Transparent
//#####################################################################
template<class T> bool OPENGL_OBJECT<T>::
Is_Transparent() const
{
    return false;  // this is usually what we want
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T> void OPENGL_OBJECT<T>::
Turn_Smooth_Shading_On()
{
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T> void OPENGL_OBJECT<T>::
Turn_Smooth_Shading_Off()
{
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_OBJECT<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    return 0;
}
//#####################################################################
// Function Set_Selection
//#####################################################################
template<class T> void OPENGL_OBJECT<T>::
Set_Selection(OPENGL_SELECTION<T>* selection)
{
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_OBJECT<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_OBJECT<T>::
Clear_Highlight()
{
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_OBJECT<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const
{
}
//#####################################################################
// Function Selection_Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_OBJECT<T>::
Selection_Bounding_Box(OPENGL_SELECTION<T>* selection) const
{
    return RANGE<TV>::Centered_Box();
}
//#####################################################################
// Function Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_OBJECT<T>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION<T>* old_selection,bool& delete_selection)
{
    return 0;
}
//#####################################################################
// Function Set_Slice
//#####################################################################
template<class T> void OPENGL_OBJECT<T>::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;Slice_Has_Changed();
}
//#####################################################################
// Function Slice_Has_Changed
//#####################################################################
template<class T> void OPENGL_OBJECT<T>::
Slice_Has_Changed()
{
}
namespace PhysBAM{
template class OPENGL_OBJECT<double>;
template class OPENGL_OBJECT<float>;
}
