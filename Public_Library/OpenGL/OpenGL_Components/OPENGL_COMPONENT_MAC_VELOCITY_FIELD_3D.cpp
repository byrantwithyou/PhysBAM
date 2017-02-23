//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Grid_Tools/Computations/VORTICITY_UNIFORM.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D(STREAM_TYPE stream_type,const GRID<TV> &grid, const std::string &velocity_filename)
    :OPENGL_COMPONENT<T>(stream_type,"MAC Velocity Field"),
    opengl_mac_velocity_field(stream_type,*new GRID<TV>(grid)),
     opengl_vorticity_magnitude(stream_type,opengl_mac_velocity_field.grid,opengl_vorticity_magnitude_array,OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1)),
    draw_vorticity(false),selected_index(-1,-1,-1),velocity_filename(velocity_filename),valid(false),min_vorticity(FLT_MAX),max_vorticity(FLT_MIN)
{
    viewer_callbacks.Set("toggle_velocity_mode",{[this](){Toggle_Velocity_Mode();},"Toggle velocity mode"});
    viewer_callbacks.Set("toggle_velocity_mode_and_draw",{[this](){Toggle_Velocity_Mode_And_Draw();},"Toggle velocity mode and draw"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_arrowhead",{[this](){Toggle_Arrowhead();},"Toggle arrowhead"});
    viewer_callbacks.Set("toggle_draw_vorticity",{[this](){Toggle_Draw_Vorticity();},"Toggle draw vorticity"});
    viewer_callbacks.Set("normalize_vorticity_color_map",{[this](){Normalize_Vorticity_Color_Map();},"Normalize vorticity map based on current frame"});

    is_animation = FILE_UTILITIES::Is_Animated(velocity_filename);
    opengl_vorticity_magnitude.Set_Scale_Range(0,100);
    frame_loaded = -1;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
~OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D()
{
    delete &opengl_mac_velocity_field.grid;
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(velocity_filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Displa
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Display() const
{
    if(valid && draw){
        opengl_mac_velocity_field.Display();
        if(draw_vorticity) opengl_vorticity_magnitude.Display();}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_mac_velocity_field.Bounding_Box();
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Reinitialize()
{
    if(draw)
    {
        if(!valid ||
            (is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(velocity_filename, frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File(stream_type,tmp_filename,opengl_mac_velocity_field.face_velocities);//u,opengl_mac_velocity_field.v,opengl_mac_velocity_field.w);
            else
                return;

            opengl_mac_velocity_field.Update();

            if(draw_vorticity){
                Update_Vorticity();
                opengl_vorticity_magnitude.Update();}

            frame_loaded = frame;
            valid = true;
        }
    }
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Print_Selection_Info(std::ostream& stream) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": ";
        opengl_mac_velocity_field.Print_Selection_Info(stream);
        if(draw_vorticity) stream<<"vorticity magnitude = "<<opengl_vorticity_magnitude.values(selected_index)<<std::endl;}
}
//#####################################################################
// Function Toggle_Velocity_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Toggle_Velocity_Mode()
{
    opengl_mac_velocity_field.Toggle_Velocity_Mode();
}
//#####################################################################
// Function Toggle_Velocity_Mode_And_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Toggle_Velocity_Mode_And_Draw()
{
    if(draw)
    {
        Toggle_Velocity_Mode();
        if((int)opengl_mac_velocity_field.velocity_mode==0) Set_Draw(!draw);
    }
    else Set_Draw(!draw);
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Increase_Vector_Size()
{
    opengl_mac_velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Decrease_Vector_Size()
{
    opengl_mac_velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Toggle_Arrowhead()
{
    opengl_mac_velocity_field.Toggle_Arrowhead_Mode();
}
//#####################################################################
// Function Toggle_Draw_Vorticity
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Toggle_Draw_Vorticity()
{
    draw_vorticity=!draw_vorticity;
    if(draw_vorticity) valid=false;
    Reinitialize();
}
//#####################################################################
// Function Normalize_Vorticity_Color_Map
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Normalize_Vorticity_Color_Map()
{
    if(!draw_vorticity) return;
    opengl_vorticity_magnitude.Set_Scale_Range(min_vorticity,max_vorticity);
    valid=false;
    Reinitialize();
}
//#####################################################################
// Function Update_Vorticity
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>::
Update_Vorticity()
{
    GRID<TV>& grid=opengl_mac_velocity_field.grid;
    RANGE<VECTOR<int,3> > domain_indices(grid.Domain_Indices());domain_indices.Change_Size(-VECTOR<int,3>::All_Ones_Vector());
    FACE_LOOKUP_UNIFORM<TV> lookup(opengl_mac_velocity_field.face_velocities);
    opengl_vorticity_magnitude.values.Resize(grid.Domain_Indices());
    for(CELL_ITERATOR<TV> iterator(grid,domain_indices);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
        T vorticity_magnitude=VORTICITY_UNIFORM<TV>::Vorticity(grid,lookup,index).Magnitude();
        opengl_vorticity_magnitude.values(index)=vorticity_magnitude;
        min_vorticity=min(min_vorticity,vorticity_magnitude);
        max_vorticity=max(max_vorticity,vorticity_magnitude);}
}
namespace PhysBAM{
template class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<float>;
template class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<double>;
}
