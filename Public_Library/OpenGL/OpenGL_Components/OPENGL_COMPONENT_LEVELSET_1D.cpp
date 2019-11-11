//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_1D
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_LEVELSET_1D
//#####################################################################
template<class T> OPENGL_COMPONENT_LEVELSET_1D<T>::
OPENGL_COMPONENT_LEVELSET_1D(const VIEWER_DIR& viewer_dir,GRID<TV> &grid,const std::string& levelset_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :OPENGL_COMPONENT<T>(viewer_dir,"Levelset 1D"),levelset_filename(levelset_filename_input),opengl_levelset(0)
{
    opengl_levelset=new OPENGL_LEVELSET_1D<T>(*(new LEVELSET<TV>(grid,*(new ARRAY<T,TV_INT>))),point_color,line_color);
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
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_LEVELSET_1D<T>::
Set_Frame()
{
    
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
    if(!draw) return;
    valid=false;
    std::string filename=viewer_dir.current_directory+"/"+levelset_filename;
    if(!File_Exists(filename)) return;
    Read_From_File(filename,opengl_levelset->levelset);
    valid=true;
}
//##################################################################### 
namespace PhysBAM{
template class OPENGL_COMPONENT_LEVELSET_1D<float>;
template class OPENGL_COMPONENT_LEVELSET_1D<double>;
}
