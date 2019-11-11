//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2>::
OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid_input,const std::string &values_filename_input,OPENGL_COLOR_MAP<T2>* color_map_input)
    :OPENGL_COMPONENT<T>(viewer_dir,"Face Scalar Field 2D"),opengl_scalar_field(grid_input,opengl_scalar_field_data,color_map_input),
    values_filename(values_filename_input),valid(false)
{
    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2>::
~OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D()
{
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2>::
Display() const
{
    if(valid && draw) opengl_scalar_field.Display();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2>::
Print_Selection_Info(std::ostream& stream) const
{
    stream<<component_name<<": "<<std::endl;
    opengl_scalar_field.Print_Selection_Info(stream);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2>::
Bounding_Box() const
{
    if(valid && draw) return opengl_scalar_field.Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T2>::
Reinitialize()
{
    if(!draw) return;
    valid = false;
    std::string filename=viewer_dir.current_directory+"/"+values_filename;
    if(File_Exists(filename))
        Read_From_File(filename,opengl_scalar_field.face_values);
    else return;
    opengl_scalar_field.Update();
    valid = true;
}
namespace PhysBAM{
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<float,int>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<float,bool>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<float,float>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<double,int>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<double,bool>;
template class OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<double,double>;
}
