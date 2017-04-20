//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T>::
OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D(STREAM_TYPE stream_type,const GRID<TV> &grid,const std::string &filename_input)
    :OPENGL_COMPONENT<T>(stream_type,"Pseudo Dirichlet"),mac_grid(grid.Get_MAC_Grid()),velocity_scale(0.025),filename(filename_input),frame_loaded(-1),valid(false)
{
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});

    is_animation = Is_Animated(filename);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T>::
Valid_Frame(int frame_input) const
{
    return Frame_File_Exists(filename, frame_input);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T>::
Display() const
{
    OPENGL_COLOR point_color=OPENGL_COLOR::White();
    OPENGL_COLOR velocity_color=OPENGL_COLOR::Yellow();
    if(valid && draw){
        glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_LINE_BIT);
        glDisable(GL_LIGHTING);glPointSize(5.0);
        T x_offset=0.4*mac_grid.dX.x,y_offset=0.4*mac_grid.dX.y;
        point_color.Send_To_GL_Pipeline();
        OpenGL_Begin(GL_POINTS);
        for(int i=0;i<pseudo_dirichlet_cells.m;i++){
            const VECTOR<int,2>& cell_index=pseudo_dirichlet_cells(i).x;char face_mask=pseudo_dirichlet_cells(i).z;
            VECTOR<T,2> pos=mac_grid.X(cell_index);
            if(face_mask&0x01) OpenGL_Vertex(VECTOR<T,2>(pos.x-x_offset,pos.y));
            if(face_mask&0x02) OpenGL_Vertex(VECTOR<T,2>(pos.x+x_offset,pos.y));
            if(face_mask&0x04) OpenGL_Vertex(VECTOR<T,2>(pos.x,pos.y-y_offset));
            if(face_mask&0x08) OpenGL_Vertex(VECTOR<T,2>(pos.x,pos.y+y_offset));}
        OpenGL_End();
        OpenGL_Begin(GL_LINES);
        velocity_color.Send_To_GL_Pipeline();
        for(int i=0;i<pseudo_dirichlet_cells.m;i++){
            const VECTOR<int,2>& cell_index=pseudo_dirichlet_cells(i).x;const VECTOR<T,2>& velocity=pseudo_dirichlet_cells(i).y;
            VECTOR<T,2> pos=mac_grid.X(cell_index);
            OPENGL_SHAPES::Draw_Arrow(pos,pos+velocity_scale*velocity);}
        OpenGL_End();
        glPopAttrib();
    }
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T>::
Reinitialize(bool force)
{
    if(draw||force)
    {
        if(!valid || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string tmp_filename = Get_Frame_Filename(filename, frame);
            if(File_Exists(tmp_filename))
                Read_From_File(stream_type,tmp_filename,pseudo_dirichlet_cells);
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
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(mac_grid.domain);
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T>::
Increase_Vector_Size()
{
    velocity_scale*=(T)1.1;
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<T>::
Decrease_Vector_Size()
{
    velocity_scale/=(T)1.1;
}
namespace PhysBAM{
template class OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<float>;
template class OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D<double>;
}
