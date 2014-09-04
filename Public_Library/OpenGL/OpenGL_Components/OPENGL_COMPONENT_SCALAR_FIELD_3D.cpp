//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
OPENGL_COMPONENT_SCALAR_FIELD_3D(const GRID<TV> &grid_input, const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input)
    :OPENGL_COMPONENT("Scalar Field 3D"), opengl_scalar_field(grid_input,*new ARRAY<T2,VECTOR<int,3> >,color_map_input),
      scalar_field_filename(scalar_field_filename_input), frame_loaded(-1), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(scalar_field_filename);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
OPENGL_COMPONENT_SCALAR_FIELD_3D(const GRID<TV> &grid_input, const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,
                                 typename OPENGL_SCALAR_FIELD_3D<T,T2>::DRAW_MODE draw_mode_input)
    :OPENGL_COMPONENT("Scalar Field 3D"), opengl_scalar_field(grid_input,*new ARRAY<T2,VECTOR<int,3> >,color_map_input,draw_mode_input),
      scalar_field_filename(scalar_field_filename_input), frame_loaded(-1), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(scalar_field_filename);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2,class RW> OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
~OPENGL_COMPONENT_SCALAR_FIELD_3D()
{
    delete &opengl_scalar_field.values;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2,class RW> bool OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(scalar_field_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    if(draw_input) opengl_scalar_field.Set_Slice(slice);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Display() const
{
    if(valid && draw) opengl_scalar_field.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_scalar_field.Print_Selection_Info(output_stream,current_selection);}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Reinitialize()
{
    if(draw)
    {
        if((is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string filename=FILE_UTILITIES::Get_Frame_Filename(scalar_field_filename,frame);
            if(FILE_UTILITIES::File_Exists(filename))
                FILE_UTILITIES::Read_From_File<RW>(filename,opengl_scalar_field.values);
            else
                return;

            opengl_scalar_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}
//#####################################################################
// Function Toggle_Smooth_Slice
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Toggle_Smooth_Slice()
{
    opengl_scalar_field.Toggle_Smooth_Slice_Texture();
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Toggle_Draw_Mode()
{
    opengl_scalar_field.Toggle_Draw_Mode();
}
//#####################################################################
// Function Toggle_Color_Map
//#####################################################################
template<class T,class T2,class RW> void OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2,RW>::
Toggle_Color_Map()
{
    opengl_scalar_field.Toggle_Color_Map();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<float,int,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<float,bool,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<float,float,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<double,int,double>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<double,bool,double>;
template class OPENGL_COMPONENT_SCALAR_FIELD_3D<double,double,double>;
}
