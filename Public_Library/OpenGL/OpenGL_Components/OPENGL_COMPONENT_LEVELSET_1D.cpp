//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_1D
//##################################################################### 
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_LEVELSET_1D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
OPENGL_COMPONENT_LEVELSET_1D(GRID<TV> &grid,const std::string& levelset_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :OPENGL_COMPONENT<T>("Levelset 1D"),levelset_filename(levelset_filename_input),opengl_levelset(0)
{
    is_animation=FILE_UTILITIES::Is_Animated(levelset_filename);
    opengl_levelset=new OPENGL_LEVELSET_1D<T>(*(new LEVELSET<TV>(grid,*(new ARRAY<T,TV_INT>))),point_color,line_color);
    Reinitialize();
}
//#####################################################################
// Function ~OPENGL_COMPONENT_LEVELSET_1D
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
~OPENGL_COMPONENT_LEVELSET_1D()
{
    delete &opengl_levelset->levelset;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(levelset_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Display() const
{
    opengl_levelset->Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return opengl_levelset->Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_LEVELSET_1D<T,RW>::
Reinitialize()
{
    if(draw && ((is_animation && frame_loaded != frame) || (!is_animation && frame_loaded<0))){
        valid=false;
        std::string filename=FILE_UTILITIES::Get_Frame_Filename(levelset_filename,frame);
        if(FILE_UTILITIES::File_Exists(filename)){
            FILE_UTILITIES::Read_From_File<RW>(filename,opengl_levelset->levelset);
            frame_loaded=frame;valid=true;}}
}
//##################################################################### 
namespace PhysBAM{
template class OPENGL_COMPONENT_LEVELSET_1D<float,float>;
template class OPENGL_COMPONENT_LEVELSET_1D<float,double>;
template class OPENGL_COMPONENT_LEVELSET_1D<double,double>;
}
