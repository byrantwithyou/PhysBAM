//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_1D
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_LEVELSET_1D
//#####################################################################
template<class T> OPENGL_COMPONENT_LEVELSET_1D<T>::
OPENGL_COMPONENT_LEVELSET_1D(STREAM_TYPE stream_type,GRID<TV> &grid,const std::string& levelset_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :OPENGL_COMPONENT<T>(stream_type,"Levelset 1D"),levelset_filename(levelset_filename_input),opengl_levelset(0)
{
    is_animation=Is_Animated(levelset_filename);
    opengl_levelset=new OPENGL_LEVELSET_1D<T>(stream_type,*(new LEVELSET<TV>(grid,*(new ARRAY<T,TV_INT>))),point_color,line_color);
    Reinitialize();
}
//#####################################################################
// Function ~OPENGL_COMPONENT_LEVELSET_1D
//#####################################################################
template<class T> OPENGL_COMPONENT_LEVELSET_1D<T>::
~OPENGL_COMPONENT_LEVELSET_1D()
{
    delete &opengl_levelset->levelset;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_LEVELSET_1D<T>::
Valid_Frame(int frame_input) const
{
    return Frame_File_Exists(levelset_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_1D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_1D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_1D<T>::
Display() const
{
    opengl_levelset->Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_LEVELSET_1D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_levelset->Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_1D<T>::
Reinitialize()
{
    if(draw && ((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))){
        valid=false;
        std::string filename=Get_Frame_Filename(levelset_filename,frame);
        if(File_Exists(filename)){
            Read_From_File(stream_type,filename,opengl_levelset->levelset);
            frame_loaded=frame;valid=true;}}
}
//##################################################################### 
namespace PhysBAM{
template class OPENGL_COMPONENT_LEVELSET_1D<float>;
template class OPENGL_COMPONENT_LEVELSET_1D<double>;
}
