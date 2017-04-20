//#####################################################################
// Copyright 2006-2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_SCALAR_FIELD_1D
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_1D.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
OPENGL_COMPONENT_SCALAR_FIELD_1D(STREAM_TYPE stream_type,const GRID<TV> &grid,const std::string &scalar_field_filename,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :OPENGL_COMPONENT<T>(stream_type,"Scalar Field 1D"),scalar_field_filename(scalar_field_filename),frame_loaded(INT_MIN),valid(false),opengl_scalar_field(stream_type,grid,*new ARRAY<T2,VECTOR<int,1> >,point_color,line_color)
{
    viewer_callbacks.Set("increase_scale",{[this](){Increase_Scale();},"Increase scale"});
    viewer_callbacks.Set("decrease_scale",{[this](){Decrease_Scale();},"Decrease scale"});
    is_animation=Is_Animated(scalar_field_filename);
}
//#####################################################################
// Function ~OPENGL_COMPONENT_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
~OPENGL_COMPONENT_SCALAR_FIELD_1D()
{
    delete &opengl_scalar_field.values;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class T2> bool OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
Valid_Frame(int frame_input) const
{
    return Frame_File_Exists(scalar_field_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
Display() const
{
    if(valid && draw) opengl_scalar_field.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
Bounding_Box() const
{
    if(valid && draw) return opengl_scalar_field.Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
Reinitialize()
{
    if(draw && ((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))){
        valid=false;
        std::string filename=Get_Frame_Filename(scalar_field_filename,frame);
        if(File_Exists(filename)){
            Read_From_File(stream_type,filename,opengl_scalar_field.values);
            frame_loaded=frame;valid=true;}}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
Scale(const T scale)
{
    opengl_scalar_field.Scale(scale);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
Increase_Scale()
{
    Scale((T)2);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T2>::
Decrease_Scale()
{
    Scale((T).5);
}
//##################################################################### 
namespace PhysBAM{
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<float,float>;
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<float,bool>;
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<double,bool>;
template class OPENGL_COMPONENT_SCALAR_FIELD_1D<double,double>;
}
