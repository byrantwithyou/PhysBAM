//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
OPENGL_COMPONENT_TRIANGULATED_SURFACE(const VIEWER_DIR& viewer_dir,const std::string& filename,bool use_display_list)
    :OPENGL_COMPONENT<T>(viewer_dir,"Triangulated Surface"), 
      triangulated_surface(*TRIANGULATED_SURFACE<T>::Create()),
    opengl_triangulated_surface(triangulated_surface, false,
                                  OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Red()),
                                  OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Blue())),
      filename(filename),  valid(false), use_display_list(use_display_list)
{
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
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
Set_Frame()
{
    
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
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
Selection_Bounding_Box() const
{
    if(valid && draw) return opengl_triangulated_surface.Selection_Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>::
Reinitialize()
{
    if(!draw) return;
    valid = false;
    std::string tmp_filename = viewer_dir.current_directory+"/"+filename;
    if(File_Exists(tmp_filename))
        Read_From_File(tmp_filename,triangulated_surface);
    else
        return;

    opengl_triangulated_surface.name=tmp_filename;
    valid = true;
}
namespace PhysBAM{
template class OPENGL_COMPONENT_TRIANGULATED_SURFACE<float>;
template class OPENGL_COMPONENT_TRIANGULATED_SURFACE<double>;
}

