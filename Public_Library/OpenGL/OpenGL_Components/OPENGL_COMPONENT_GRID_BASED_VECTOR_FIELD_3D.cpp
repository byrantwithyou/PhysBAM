//#####################################################################
// Copyright 2004, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid,const std::string &vector_field_filename)
    :OPENGL_COMPONENT<T>(viewer_dir,"Cell Centered Velocity Field"), 
    opengl_grid_based_vector_field(*(new GRID<TV>(grid)), *(new ARRAY<VECTOR<T,3> ,VECTOR<int,3> >)), 
      vector_field_filename(vector_field_filename), valid(false)
{
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_arrowhead",{[this](){Toggle_Arrowhead();},"Toggle arrowhead"});
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
~OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D()
{
    delete &opengl_grid_based_vector_field.grid;
    delete &opengl_grid_based_vector_field.V;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
Display() const
{
    if(valid && draw) opengl_grid_based_vector_field.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(opengl_grid_based_vector_field.Bounding_Box());
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
Reinitialize()
{
    if(!draw) return;
    valid = false;
    std::string tmp_filename = viewer_dir.current_directory+"/"+vector_field_filename;
    if(File_Exists(tmp_filename))
        Read_From_File(tmp_filename,opengl_grid_based_vector_field.V);
    else
        return;
    opengl_grid_based_vector_field.Update();
    valid = true;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
Print_Selection_Info(std::ostream& stream) const
{
    stream<<component_name<<": ";
    opengl_grid_based_vector_field.Print_Selection_Info(stream);
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
Increase_Vector_Size()
{
    opengl_grid_based_vector_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
Decrease_Vector_Size()
{
    opengl_grid_based_vector_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>::
Toggle_Arrowhead()
{
    opengl_grid_based_vector_field.Toggle_Arrowhead_Mode();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<float>;
template class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<double>;
}
