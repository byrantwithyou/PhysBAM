//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/MATRIX_2X2.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Reinitialize
//##################################################################### 
template<class T> void OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<T>::
Reinitialize(bool force_load_even_if_not_drawn)
{
    if(!draw&&!force_load_even_if_not_drawn)return;
    if(!valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
        valid=false;std::string tmp_filename=FILE_UTILITIES::Get_Frame_Filename(field_filename,frame);
        if(FILE_UTILITIES::File_Exists(tmp_filename))FILE_UTILITIES::Read_From_File(stream_type,tmp_filename,field);else return;
        opengl_symmetric_matrix_field.Update();
        frame_loaded=frame;valid=true;}
}
//##################################################################### 
namespace PhysBAM{
template class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<float>;
template class OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<double>;
}
