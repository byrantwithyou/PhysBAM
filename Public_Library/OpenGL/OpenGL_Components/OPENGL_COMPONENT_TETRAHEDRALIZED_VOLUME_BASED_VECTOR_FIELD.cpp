//#####################################################################
// Copyright 2004, Zhaosheng Bao.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD(STREAM_TYPE stream_type,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const std::string &vector_field_filename)
    :OPENGL_COMPONENT<T>(stream_type,"Triangulated Area Based Vector Field 3D"),opengl_vector_field(stream_type,tetrahedralized_volume,*new ARRAY<VECTOR<T,3> >), 
    vector_field_filename(vector_field_filename),frame_loaded(-1),valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(vector_field_filename);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
~OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD()
{
    delete &opengl_vector_field.V;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(vector_field_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Display() const
{
    if(valid && draw) opengl_vector_field.Display();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_vector_field.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
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
                FILE_UTILITIES::Read_From_File(stream_type,tmp_filename,opengl_vector_field.V);
            else
                return;

            opengl_vector_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Increase_Vector_Size()
{
    opengl_vector_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Decrease_Vector_Size()
{
    opengl_vector_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>::
Toggle_Arrowhead()
{
    opengl_vector_field.draw_arrowhead = !opengl_vector_field.draw_arrowhead;
}
namespace PhysBAM{
template class OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<float>;
template class OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<double>;
}
