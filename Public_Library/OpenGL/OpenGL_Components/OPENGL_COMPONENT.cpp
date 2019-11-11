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
OPENGL_COMPONENT(const VIEWER_DIR& viewer_dir,const std::string &name)
    :viewer_dir(viewer_dir),draw(true),component_name(name)
{
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
