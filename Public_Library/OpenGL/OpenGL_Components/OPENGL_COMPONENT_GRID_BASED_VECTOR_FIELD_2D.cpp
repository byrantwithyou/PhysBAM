//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D(STREAM_TYPE stream_type,const GRID<TV> &grid,const std::string &vector_field_filename)
    :OPENGL_COMPONENT<T>(stream_type,"Grid Based Vector Field 2D"), opengl_grid_based_vector_field(stream_type,*(new GRID<TV>(grid)), *(new ARRAY<VECTOR<T,2> ,VECTOR<int,2> >)), 
      vector_field_filename(vector_field_filename), valid(false)
{
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_arrowhead",{[this](){Toggle_Arrowhead();},"Toggle arrowhead"});

    is_animation = FILE_UTILITIES::Is_Animated(vector_field_filename);
    frame_loaded = -1;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
~OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D()
{
    delete &opengl_grid_based_vector_field.grid;
    delete &opengl_grid_based_vector_field.V;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(vector_field_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Display() const
{
    if(valid && draw) opengl_grid_based_vector_field.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_grid_based_vector_field.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": ";
        opengl_grid_based_vector_field.Print_Selection_Info(stream,selection);}
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Reinitialize(bool force_load_even_if_not_drawn)
{
    if(draw||force_load_even_if_not_drawn)
    {
        if(!valid ||
            (is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(vector_field_filename, frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File(stream_type,tmp_filename,opengl_grid_based_vector_field.V);
            else
                return;

            opengl_grid_based_vector_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Increase_Vector_Size()
{
    opengl_grid_based_vector_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Decrease_Vector_Size()
{
    opengl_grid_based_vector_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>::
Toggle_Arrowhead()
{
    opengl_grid_based_vector_field.draw_arrowhead = !opengl_grid_based_vector_field.draw_arrowhead;
}
namespace PhysBAM{
template class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<float>;
template class OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<double>;
}
