//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D(STREAM_TYPE stream_type,OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>& V_minus_component,OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>& V_plus_component,
       OPENGL_COMPONENT_LEVELSET_2D<T>& levelset_component)
    :OPENGL_COMPONENT<T>(stream_type,"Two-Phase Magnitude Velocity Field 2D"),magnitude_height_scale(0),V_minus_component(V_minus_component),V_plus_component(V_plus_component),levelset_component(levelset_component),
    opengl_two_phase_velocity_magnitude(stream_type,V_minus_component.opengl_grid_based_vector_field.grid,V_minus_component.opengl_grid_based_vector_field.V,V_plus_component.opengl_grid_based_vector_field.V,levelset_component.opengl_levelset->levelset)
{
    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
~OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D()
{}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Valid_Frame(int frame_input) const
{
    return valid;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Display() const
{
    if(valid && draw)
        opengl_two_phase_velocity_magnitude.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_two_phase_velocity_magnitude.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Toggle_3D_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Toggle_3D_Mode()
{
    if(magnitude_height_scale>0)magnitude_height_scale=0;
    else magnitude_height_scale=(T)1;
    opengl_two_phase_velocity_magnitude.Scale_Height(magnitude_height_scale);
    Reinitialize();
}
//#####################################################################
// Function Increase_Point_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Increase_Point_Size()
{
    opengl_two_phase_velocity_magnitude.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Point_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Decrease_Point_Size()
{
    opengl_two_phase_velocity_magnitude.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Reinitialize(const bool force_even_if_not_drawn)
{  
    if(draw||force_even_if_not_drawn){
        V_minus_component.Reinitialize(true);
        V_plus_component.Reinitialize(true);
        levelset_component.Reinitialize(true);
        valid=V_minus_component.Valid_Frame(frame)&&V_plus_component.Valid_Frame(frame)&&levelset_component.Valid_Frame(frame);
        opengl_two_phase_velocity_magnitude.Update();
    }
}
namespace PhysBAM{
template class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<float>;
template class OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<double>;
}
