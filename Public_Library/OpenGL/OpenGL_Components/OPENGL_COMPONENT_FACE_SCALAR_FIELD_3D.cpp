//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>::
OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D(STREAM_TYPE stream_type,const GRID<TV> &grid_input,const std::string &values_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input)
    :OPENGL_COMPONENT<T>(stream_type,"Face Scalar Field 3D"),opengl_scalar_field(stream_type,grid_input,internal_scalar_field,color_map_input),
      values_filename(values_filename_input),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(values_filename);
    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>::
~OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D()
{
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2> bool OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(values_filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    if(draw_input) opengl_scalar_field.Set_Slice(slice);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>::
Display() const
{
    if(valid && draw) opengl_scalar_field.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>::
Bounding_Box() const
{
    if(valid && draw) return opengl_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>::
Print_Selection_Info(std::ostream& stream) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        opengl_scalar_field.Print_Selection_Info(stream);}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>::
Reinitialize()
{
    if(draw)
    {
        if((is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string filename;
            filename=FILE_UTILITIES::Get_Frame_Filename(values_filename,frame);
            if(FILE_UTILITIES::File_Exists(filename))
                FILE_UTILITIES::Read_From_File(stream_type,filename,opengl_scalar_field.face_values);
            else return;

            opengl_scalar_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}
namespace PhysBAM{
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<float,int>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<float,bool>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<float,float>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<double,int>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<double,bool>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<double,double>;
}
