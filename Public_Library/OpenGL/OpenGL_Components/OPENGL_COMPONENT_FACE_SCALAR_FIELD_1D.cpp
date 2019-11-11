//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid_input,const std::string &values_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color)
    :OPENGL_COMPONENT<T>(viewer_dir,"Face Scalar Field 1D"),opengl_face_scalar_field(grid_input,*new ARRAY<T2,FACE_INDEX<1> >,point_color,line_color),
      values_filename(values_filename_input),valid(false)
{
    viewer_callbacks.Set("increase_scale",{[this](){Increase_Scale();},"Increase scale"});
    viewer_callbacks.Set("decrease_scale",{[this](){Decrease_Scale();},"Decrease scale"});
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
// Function Set_Frame
//#####################################################################
template<class T,class T2> void OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Set_Frame()
{
    
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
Print_Selection_Info(std::ostream& stream) const
{
    stream<<component_name<<": "<<std::endl;
    opengl_face_scalar_field.Print_Selection_Info(stream);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T2>::
Bounding_Box() const
{
    if(valid && draw) return opengl_face_scalar_field.Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
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
    if(!draw) return;
    valid=false;
    std::string filename=viewer_dir.current_directory+"/"+values_filename;
    if(File_Exists(filename))
        Read_From_File(filename,opengl_face_scalar_field.face_values);
    else return;
    valid = true;
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
