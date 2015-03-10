//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D(STREAM_TYPE stream_type,const GRID<TV> &grid_input,const std::string &values_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :OPENGL_COMPONENT<T>(stream_type,"Face Scalar Field 1D"),opengl_face_scalar_field(stream_type,grid_input,*new ARRAY<T2,FACE_INDEX<1> >,point_color,line_color),
      values_filename(values_filename_input),frame_loaded(-1),valid(false)
{
    viewer_callbacks.Set("increase_scale",{[this](){Increase_Scale();},"Increase scale"});
    viewer_callbacks.Set("decrease_scale",{[this](){Decrease_Scale();},"Decrease scale"});
    is_animation = FILE_UTILITIES::Is_Animated(values_filename);
    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
~OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D()
{
    delete &opengl_face_scalar_field.face_values;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2> bool OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(values_filename,frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Display() const
{
    if(valid && draw) opengl_face_scalar_field.Display();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        opengl_face_scalar_field.Print_Selection_Info(stream,selection);}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Bounding_Box() const
{
    if(valid && draw) return opengl_face_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Scale
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Scale(const T scale)
{
    opengl_face_scalar_field.Scale(scale);
}
//#####################################################################
// Function Increase_Scale
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Increase_Scale()
{
    Scale((T)2);
}
//#####################################################################
// Function Decrease_Scale 
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Decrease_Scale()
{
    Scale((T).5);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Reinitialize()
{
    if(draw && ((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))){
        valid=false;

        std::string filename=FILE_UTILITIES::Get_Frame_Filename(values_filename,frame);
        if(FILE_UTILITIES::File_Exists(filename))
            FILE_UTILITIES::Read_From_File(stream_type,filename,opengl_face_scalar_field.face_values);
        else return;
        frame_loaded = frame;
        valid = true;}
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<float,int>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<float,bool>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<float,float>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<double,int>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<double,bool>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<double,double>;
}
