//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT<T>::
OPENGL_COMPONENT(STREAM_TYPE stream_type,const std::string &name)
 :OPENGL_OBJECT<T>(stream_type),frame(0),draw(true),is_animation(false),component_name(name)
{
    viewer_callbacks.Set("next_frame",{[this](){Set_Frame(frame+1);},"Next frame"});
    viewer_callbacks.Set("prev_frame",{[this](){Set_Frame(frame-1);},"Prev frame"});
    viewer_callbacks.Set("toggle_draw",{[this](){Set_Draw(!draw);},"Toggle draw"});
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT<T>::
~OPENGL_COMPONENT()
{
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT<T>::
Valid_Frame(int frame_input) const
{
    return false;
}
//#####################################################################
// Function Is_Up_To_Date
//#####################################################################
template<class T> bool OPENGL_COMPONENT<T>::
Is_Up_To_Date(int frame) const
{
    return true;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT<T>::
Set_Frame(int frame_input)
{
    frame = frame_input;
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT<T>::
Set_Draw(bool draw_input)
{
    draw = draw_input;
}
//#####################################################################
// Function Draw_All_Objects
//#####################################################################
template<class T> void OPENGL_COMPONENT<T>::
Draw_All_Objects()
{
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT<T>::
Use_Bounding_Box() const
{
    return draw;
}
namespace PhysBAM{
template class OPENGL_COMPONENT<double>;
template class OPENGL_COMPONENT<float>;
}
