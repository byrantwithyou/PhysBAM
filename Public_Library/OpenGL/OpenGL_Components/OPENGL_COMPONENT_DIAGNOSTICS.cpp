//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DIAGNOSTICS.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// OPENGL_COMPONENT_DIAGNOSTICS
//#####################################################################
template<class T> OPENGL_COMPONENT_DIAGNOSTICS<T>::
OPENGL_COMPONENT_DIAGNOSTICS(const VIEWER_DIR& viewer_dir,const std::string& filename_input)
    :OPENGL_COMPONENT<T>(viewer_dir),filename(filename_input),valid(false)
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
// Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_DIAGNOSTICS<T>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_DIAGNOSTICS<T>::
Print_Selection_Info(std::ostream& output_stream) const
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
    if(!draw) return;
    lines.Remove_All();
    valid=false;
    std::string tmp_filename = viewer_dir.current_directory+"/"+filename;
    if(!File_Exists(tmp_filename)) return;
    std::istream* input=Safe_Open_Input_Raw(tmp_filename,false);
    std::string line;
    while(std::getline(*input,line)) lines.Append(line);
    delete input;
    valid=true;
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
