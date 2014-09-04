//#####################################################################
// Copyright 2004, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D(const GRID<TV> &grid,const std::string &vector_field_filename)
    :OPENGL_COMPONENT<T>("Cell Centered Velocity Field"), 
      opengl_grid_based_vector_field(*(new GRID<TV>(grid)), *(new ARRAY<VECTOR<T,3> ,VECTOR<int,3> >)), 
      vector_field_filename(vector_field_filename), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(vector_field_filename);
    frame_loaded = -1;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
~OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D()
{
    delete &opengl_grid_based_vector_field.grid;
    delete &opengl_grid_based_vector_field.V;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(vector_field_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Display() const
{
    if(valid && draw) opengl_grid_based_vector_field.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(opengl_grid_based_vector_field.Bounding_Box());
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Reinitialize()
{
    if(draw)
    {
        if(!valid ||
            (is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(vector_field_filename, frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_grid_based_vector_field.V);
            else
                return;

            opengl_grid_based_vector_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": ";
        opengl_grid_based_vector_field.Print_Selection_Info(stream,selection);}
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Set_Vector_Size(double size)
{
    opengl_grid_based_vector_field.size = size;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Increase_Vector_Size()
{
    opengl_grid_based_vector_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Decrease_Vector_Size()
{
    opengl_grid_based_vector_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T,RW>::
Toggle_Arrowhead()
{
    opengl_grid_based_vector_field.Toggle_Arrowhead_Mode();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<float,float>;
template class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<double,double>;
}
