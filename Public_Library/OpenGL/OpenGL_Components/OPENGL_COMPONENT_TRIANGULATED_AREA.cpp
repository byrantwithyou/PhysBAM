//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_AREA.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TRIANGULATED_AREA<T>::
OPENGL_COMPONENT_TRIANGULATED_AREA(const VIEWER_DIR& viewer_dir,const std::string &filename)
    :OPENGL_COMPONENT<T>(viewer_dir,"Triangulated Surface"),color_map(0), 
      triangulated_area(*TRIANGULATED_AREA<T>::Create()),
    opengl_triangulated_area(triangulated_area),
      filename(filename),color_map_filename(0),valid(false)
{
    Reinitialize();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TRIANGULATED_AREA<T>::
OPENGL_COMPONENT_TRIANGULATED_AREA(const VIEWER_DIR& viewer_dir,const std::string &filename,const std::string &color_map_filename_input)
    :OPENGL_COMPONENT<T>(viewer_dir,"Triangulated Surface"),color_map(new ARRAY<OPENGL_COLOR >), 
      triangulated_area(*TRIANGULATED_AREA<T>::Create()),
    opengl_triangulated_area(triangulated_area,false,OPENGL_COLOR::Red(),OPENGL_COLOR::Black()),
      filename(filename),color_map_filename(&color_map_filename_input),valid(false)
{
    opengl_triangulated_area.Set_Color_Map(color_map);
    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TRIANGULATED_AREA<T>::
~OPENGL_COMPONENT_TRIANGULATED_AREA()
{
    delete &triangulated_area.mesh;
    delete &triangulated_area.particles;
    delete &triangulated_area;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_AREA<T>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_AREA<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_AREA<T>::
Display() const
{
    if(valid && draw) opengl_triangulated_area.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_TRIANGULATED_AREA<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_triangulated_area.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_TRIANGULATED_AREA<T>::
Reinitialize()
{
    if(!draw) return;
    valid = false;
    std::string tmp_filename = viewer_dir.current_directory+"/"+filename;
    if(File_Exists(tmp_filename))
        Read_From_File(tmp_filename,triangulated_area);
    else
        return;
    if(color_map) {
        std::string tmp_color_map_filename = viewer_dir.current_directory+"/"+*color_map_filename;
        //if(File_Exists(tmp_filename))
        if(File_Exists(tmp_color_map_filename))
            Read_From_File(tmp_color_map_filename,*color_map);
        else
            return;
    }            
    valid = true;
}
namespace PhysBAM{
template class OPENGL_COMPONENT_TRIANGULATED_AREA<float>;
template class OPENGL_COMPONENT_TRIANGULATED_AREA<double>;
}
