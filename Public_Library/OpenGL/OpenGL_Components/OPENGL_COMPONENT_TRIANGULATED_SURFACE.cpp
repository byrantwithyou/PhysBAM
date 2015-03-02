//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
OPENGL_COMPONENT_TRIANGULATED_SURFACE(STREAM_TYPE stream_type,const std::string &filename, bool use_display_list)
    :OPENGL_COMPONENT<T>(stream_type,"Triangulated Surface"), 
      triangulated_surface(*TRIANGULATED_SURFACE<T>::Create()),
    opengl_triangulated_surface(stream_type,triangulated_surface, false,
                                  OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Red()),
                                  OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Blue())),
      filename(filename), frame_loaded(-1), valid(false), use_display_list(use_display_list)
{
    is_animation = FILE_UTILITIES::Is_Animated(filename);
    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
~OPENGL_COMPONENT_TRIANGULATED_SURFACE()
{
    delete &triangulated_surface.mesh;
    delete &triangulated_surface.particles;
    delete &triangulated_surface;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
Display() const
{
    if(valid && draw) 
    {
        if(slice && slice->Is_Slice_Mode()) {
            glPushAttrib(GL_ENABLE_BIT);
            slice->Enable_Clip_Planes();
        }
        opengl_triangulated_surface.Display();
        if(slice && slice->Is_Slice_Mode()) {
            glPopAttrib();
        }
    }
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_triangulated_surface.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
Reinitialize()
{
    if(draw)
    {
        if((is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;
            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(filename, frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File(stream_type,tmp_filename,triangulated_surface);
            else
                return;

            opengl_triangulated_surface.name=tmp_filename;
            frame_loaded = frame;
            valid = true;
        }
    }
}
namespace PhysBAM{
template class OPENGL_COMPONENT_TRIANGULATED_SURFACE<float>;
template class OPENGL_COMPONENT_TRIANGULATED_SURFACE<double>;
}

