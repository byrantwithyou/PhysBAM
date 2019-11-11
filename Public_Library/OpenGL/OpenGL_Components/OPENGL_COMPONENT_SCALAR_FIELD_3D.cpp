//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
OPENGL_COMPONENT_SCALAR_FIELD_3D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid_input,
    const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,
    typename OPENGL_SCALAR_FIELD_3D<T,T2>::DRAW_MODE draw_mode_input)
    :OPENGL_COMPONENT<T>(viewer_dir,"Scalar Field 3D"), opengl_scalar_field(grid_input,*new ARRAY<T2,VECTOR<int,3> >,color_map_input,draw_mode_input),
      scalar_field_filename(scalar_field_filename_input),  valid(false)
{
    viewer_callbacks.Set("toggle_smooth_slice",{[this](){Toggle_Smooth_Slice();},"Toggle smooth"});
    viewer_callbacks.Set("toggle_draw_mode",{[this](){Toggle_Draw_Mode();},"Toggle draw mode"});
    viewer_callbacks.Set("toggle_color_map",{[this](){Toggle_Color_Map();},"Toggle color map"});
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
~OPENGL_COMPONENT_SCALAR_FIELD_3D()
{
    delete &opengl_scalar_field.values;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    if(draw_input) opengl_scalar_field.Set_Slice(slice);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
Display() const
{
    if(valid && draw) opengl_scalar_field.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
Bounding_Box() const
{
    if(valid && draw) return opengl_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
Print_Selection_Info(std::ostream& output_stream) const
{
    if(valid && draw){
        output_stream<<component_name<<": ";
        opengl_scalar_field.Print_Selection_Info(output_stream);}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
Reinitialize()
{
    if(!draw) return;
    valid=false;
    std::string filename=viewer_dir.current_directory+"/"+scalar_field_filename;
    if(!File_Exists(filename)) return;
    Read_From_File(filename,opengl_scalar_field.values);
    opengl_scalar_field.Update();
    valid=true;
}
//#####################################################################
// Function Toggle_Smooth_Slice
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
Toggle_Smooth_Slice()
{
    opengl_scalar_field.Toggle_Smooth_Slice_Texture();
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
Toggle_Draw_Mode()
{
    opengl_scalar_field.Toggle_Draw_Mode();
}
//#####################################################################
// Function Toggle_Color_Map
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>::
Toggle_Color_Map()
{
    opengl_scalar_field.Toggle_Color_Map();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<float,int>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<float,bool>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<float,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<double,int>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<double,bool>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<double,double>;
}
