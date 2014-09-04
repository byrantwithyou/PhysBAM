//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DIAGNOSTICS.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// OPENGL_COMPONENT_DIAGNOSTICS
//#####################################################################
template<class T> OPENGL_COMPONENT_DIAGNOSTICS<T>::
OPENGL_COMPONENT_DIAGNOSTICS(const std::string& filename_input)
    :filename(filename_input),frame_loaded(INT_MIN),valid(false)
{
}
//#####################################################################
// ~OPENGL_COMPONENT_DIAGNOSTICS
//#####################################################################
template<class T> OPENGL_COMPONENT_DIAGNOSTICS<T>::
~OPENGL_COMPONENT_DIAGNOSTICS()
{
}
//#####################################################################
// Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DIAGNOSTICS<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename,frame_input);
}
//#####################################################################
// Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_DIAGNOSTICS<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_DIAGNOSTICS<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* selection) const
{
    output_stream<<component_name<<":"<<std::endl;
    for(int i=0;i<lines.m;i++) output_stream<<"   "<<lines(i)<<std::endl;
}
//#####################################################################
// Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DIAGNOSTICS<T>::
Reinitialize()
{
    if(draw){
      if((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded==INT_MIN)){
            lines.Remove_All();
            valid=false;
            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(filename, frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename)){
                std::istream* input=FILE_UTILITIES::Safe_Open_Input(tmp_filename,false);
                std::string line;while(std::getline(*input,line)) lines.Append(line);
                delete input;frame_loaded=frame;valid=true;}}}
}
//#####################################################################
// Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DIAGNOSTICS<T>::
Use_Bounding_Box() const
{
    return false;
}
//#####################################################################
// Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_DIAGNOSTICS<T>::
Display() const
{
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_DIAGNOSTICS<double>;
template class OPENGL_COMPONENT_DIAGNOSTICS<float>;
}
