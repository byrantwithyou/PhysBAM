#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD.h>
using namespace PhysBAM;

template<class T_GRID,class RW> OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD(T_GRID& grid,const std::string& scalar_field_filename_input,OPENGL_COLOR_MAP<T>* color_map_input)
: OPENGL_COMPONENT("Adaptive Scalar Field"),opengl_adaptive_node_scalar_field(grid,*(new ARRAY<T>),color_map_input), 
                    scalar_field_filename(scalar_field_filename_input), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(scalar_field_filename);
    frame_loaded = -1;

    Reinitialize();
}

template<class T_GRID,class RW> OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
~OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD()
{
}

template<class T_GRID,class RW> bool OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(scalar_field_filename, frame_input);
}

template<class T_GRID,class RW> void OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T_GRID,class RW> void OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T_GRID,class RW> void OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_adaptive_node_scalar_field.Display(in_color);
}

template<class T_GRID,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_adaptive_node_scalar_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T_GRID,class RW> void OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Reinitialize()
{
    if (draw)
    {
        if (!valid ||
            (is_animation && (frame_loaded != frame)) ||
            (!is_animation && (frame_loaded < 0)))
        {
            valid = false;

            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(scalar_field_filename, frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_adaptive_node_scalar_field.value);
            else
                return;

            opengl_adaptive_node_scalar_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}

template<class T_GRID,class RW> void OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_adaptive_node_scalar_field.Print_Selection_Info(output_stream,current_selection);}
}

template<class T_GRID,class RW> void OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Increase_Point_Size()
{
    opengl_adaptive_node_scalar_field.point_size++;
}

template<class T_GRID,class RW> void OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Decrease_Point_Size()
{
    opengl_adaptive_node_scalar_field.point_size--;
}

template<class T_GRID,class RW> void OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID,RW>::
Toggle_Smooth()
{
    opengl_adaptive_node_scalar_field.Set_Smooth_Shading(!opengl_adaptive_node_scalar_field.smooth_shading);
}

template class OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<float>,float>;
template class OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<RED_GREEN_GRID_2D<float>,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<double>,double>;
template class OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<RED_GREEN_GRID_2D<double>,double>;
#endif
#endif
