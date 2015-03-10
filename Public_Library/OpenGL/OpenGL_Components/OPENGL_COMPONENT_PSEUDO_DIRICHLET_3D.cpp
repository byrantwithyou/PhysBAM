//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D(STREAM_TYPE stream_type,const GRID<TV> &grid,const std::string &filename_input)
    :OPENGL_COMPONENT<T>(stream_type,"Pseudo Dirichlet"),mac_grid(grid.Get_MAC_Grid()),velocity_scale(0.025),filename(filename_input),frame_loaded(-1),valid(false)
{
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});

    is_animation = FILE_UTILITIES::Is_Animated(filename);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Display() const
{
    OPENGL_COLOR point_color=OPENGL_COLOR::White();
    OPENGL_COLOR velocity_color=OPENGL_COLOR::Yellow();
    if(valid && draw){
        glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_LINE_BIT);
        if(slice && slice->Is_Slice_Mode()) slice->Enable_Clip_Planes();
        glDisable(GL_LIGHTING);glPointSize(5.0);
        T x_offset=0.4*mac_grid.dX.x,y_offset=0.4*mac_grid.dX.y,z_offset=0.4*mac_grid.dX.z;
        point_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_POINTS);
        for(int i=0;i<pseudo_dirichlet_cells.m;i++){
            const VECTOR<int,3>& cell_index=pseudo_dirichlet_cells(i).x;char face_mask=pseudo_dirichlet_cells(i).z;
            VECTOR<T,3> pos=mac_grid.X(cell_index);
            if(face_mask&0x01) OpenGL_Vertex(VECTOR<T,3>(pos.x-x_offset,pos.y,pos.z));
            if(face_mask&0x02) OpenGL_Vertex(VECTOR<T,3>(pos.x+x_offset,pos.y,pos.z));
            if(face_mask&0x04) OpenGL_Vertex(VECTOR<T,3>(pos.x,pos.y-y_offset,pos.z));
            if(face_mask&0x08) OpenGL_Vertex(VECTOR<T,3>(pos.x,pos.y+y_offset,pos.z));
            if(face_mask&0x10) OpenGL_Vertex(VECTOR<T,3>(pos.x,pos.y,pos.z-z_offset));
            if(face_mask&0x20) OpenGL_Vertex(VECTOR<T,3>(pos.x,pos.y,pos.z+z_offset));}
        OpenGL_End();
        velocity_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int i=0;i<pseudo_dirichlet_cells.m;i++){
            const VECTOR<int,3>& cell_index=pseudo_dirichlet_cells(i).x;const VECTOR<T,3>& velocity=pseudo_dirichlet_cells(i).y;
            VECTOR<T,3> pos=mac_grid.X(cell_index);
            OpenGL_Line(pos,pos+velocity_scale*velocity);}
        OpenGL_End();
        glPopAttrib();
    }
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Reinitialize(bool force)
{
    if(draw||force)
    {
        if(!valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(filename, frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File(stream_type,tmp_filename,pseudo_dirichlet_cells);
            else
                return;

            frame_loaded=frame;
            valid=true;
        }
    }
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Bounding_Box() const
{
    if(valid && draw) return RANGE<VECTOR<T,3> >(mac_grid.domain);
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION<T>::GRID_CELL_3D && Is_Up_To_Date(frame)){
        VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)current_selection)->index;
        // TODO: This is not an efficient lookup...
        for(int i=0;i<pseudo_dirichlet_cells.m;i++){
            if(pseudo_dirichlet_cells(i).x==index){
                output_stream<<component_name<<":  velocity = "<<pseudo_dirichlet_cells(i).y<<std::endl;}}}
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Set_Vector_Size(const T vector_size)
{
    velocity_scale=vector_size;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Increase_Vector_Size()
{
    velocity_scale*=(T)1.1;
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>::
Decrease_Vector_Size()
{
    velocity_scale/=(T)1.1;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<float>;
template class OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<double>;
}
