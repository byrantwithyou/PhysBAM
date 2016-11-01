//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
OPENGL_COMPONENT_SCALAR_FIELD_2D(STREAM_TYPE stream_type,GRID<TV> &grid_input, const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input)
    :OPENGL_COMPONENT<T>(stream_type,"Scalar Field 2D"), opengl_scalar_field(stream_type,grid_input,*new ARRAY<T2,VECTOR<int,2> >,color_map_input),
      scalar_field_filename(scalar_field_filename_input), frame_loaded(-1), valid(false)
{
    viewer_callbacks.Set("toggle_smooth",{[this](){Toggle_Smooth();},"Toggle smooth"});
    viewer_callbacks.Set("toggle_draw_mode",{[this](){Toggle_Draw_Mode();},"Toggle draw mode"});
    viewer_callbacks.Set("toggle_color_map",{[this](){Toggle_Color_Map();},"Toggle color map"});

    is_animation = FILE_UTILITIES::Is_Animated(scalar_field_filename);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
OPENGL_COMPONENT_SCALAR_FIELD_2D(STREAM_TYPE stream_type,GRID<TV> &grid_input, const std::string &scalar_field_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input,
                                 typename OPENGL_SCALAR_FIELD_2D<T,T2>::DRAW_MODE draw_mode_input)
    :OPENGL_COMPONENT<T>(stream_type,"Scalar Field 2D"), opengl_scalar_field(stream_type,grid_input,*new ARRAY<T2,VECTOR<int,2> >,color_map_input,draw_mode_input),
      scalar_field_filename(scalar_field_filename_input), frame_loaded(-1), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(scalar_field_filename);
}

//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
~OPENGL_COMPONENT_SCALAR_FIELD_2D()
{
    delete &opengl_scalar_field.values;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2> bool OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(scalar_field_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Display() const
{
    if(valid && draw) opengl_scalar_field.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Bounding_Box() const
{
    if(valid && draw) return opengl_scalar_field.Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Print_Selection_Info(std::ostream& output_stream) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_scalar_field.Print_Selection_Info(output_stream);}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Reinitialize()
{
    if(draw)
    {
        if((is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;
            std::string filename = FILE_UTILITIES::Get_Frame_Filename(scalar_field_filename, frame);
            if(FILE_UTILITIES::File_Exists(filename))
                FILE_UTILITIES::Read_From_File(stream_type,filename,opengl_scalar_field.values);
            else
                return;

            opengl_scalar_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}
//#####################################################################
// Function Toggle_Smooth
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Toggle_Smooth()
{
    opengl_scalar_field.Toggle_Smooth_Texture();
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Toggle_Draw_Mode()
{
    opengl_scalar_field.Toggle_Draw_Mode();
}
//#####################################################################
// Function Toggle_Color_Map
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_2D<T,T2>::
Toggle_Color_Map()
{
    opengl_scalar_field.Toggle_Color_Map();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<float,int>;
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<float,bool>;
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<float,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<double,int>;
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<double,bool>;
template class OPENGL_COMPONENT_SCALAR_FIELD_2D<double,double>;
}
