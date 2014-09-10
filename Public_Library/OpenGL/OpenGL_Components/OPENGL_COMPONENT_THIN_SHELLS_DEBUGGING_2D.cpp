//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D(GRID<TV> &grid,const std::string& directory)
    :OPENGL_COMPONENT<T>("Thin Shells Debugging"),grid(grid),
    invalid_color_map(OPENGL_COLOR::Red()),
    opengl_density_valid_mask(grid,density_valid_mask,&invalid_color_map,OPENGL_SCALAR_FIELD_2D<T,bool>::DRAW_POINTS),
    opengl_phi_valid_mask(grid,phi_valid_mask,&invalid_color_map,OPENGL_SCALAR_FIELD_2D<T,bool>::DRAW_POINTS),
    directory(directory),frame_loaded(-1),valid(false),
    draw_grid_visibility(false),draw_density_valid_mask(false),draw_phi_valid_mask(false)
{
    is_animation=true;
    mac_grid=grid.Get_MAC_Grid();
    u_grid=grid.Get_Face_Grid(0);
    v_grid=grid.Get_Face_Grid(1);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    // TODO: make more accurate
    return false;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
Display() const
{
    OPENGL_COLOR node_neighbor_not_visible_color=OPENGL_COLOR::Magenta(0.5,0.5);
    OPENGL_COLOR face_corners_not_visible_from_face_center_color=OPENGL_COLOR::Magenta(1);
    if(valid && draw){
        if(draw_grid_visibility){
            glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT | GL_CURRENT_BIT);
            glDisable(GL_LIGHTING);

            face_corners_not_visible_from_face_center_color.Send_To_GL_Pipeline();
            glLineWidth(1);
            OpenGL_Begin(GL_LINES);
            for(int i=face_corners_visible_from_face_center_u.domain.min_corner.x;i<face_corners_visible_from_face_center_u.domain.max_corner.x;i++) for(int j=face_corners_visible_from_face_center_u.domain.min_corner.y;j<face_corners_visible_from_face_center_u.domain.max_corner.y;j++){
                if(!face_corners_visible_from_face_center_u(i,j)(0).x){OpenGL_Line(u_grid.X(TV_INT(i,j)),grid.X(TV_INT(i,j)));}
                if(!face_corners_visible_from_face_center_u(i,j)(1).x){OpenGL_Line(u_grid.X(TV_INT(i,j)),grid.X(TV_INT(i,j+1)));}}
            for(int i=face_corners_visible_from_face_center_v.domain.min_corner.x;i<face_corners_visible_from_face_center_v.domain.max_corner.x;i++) for(int j=face_corners_visible_from_face_center_v.domain.min_corner.y;j<face_corners_visible_from_face_center_v.domain.max_corner.y;j++){
                if(!face_corners_visible_from_face_center_v(i,j)(0).x){OpenGL_Line(v_grid.X(TV_INT(i,j)),grid.X(TV_INT(i,j)));}
                if(!face_corners_visible_from_face_center_v(i,j)(1).x){OpenGL_Line(v_grid.X(TV_INT(i,j)),grid.X(TV_INT(i+1,j)));}}
            OpenGL_End();

            node_neighbor_not_visible_color.Send_To_GL_Pipeline();
            glLineWidth(5);
            OpenGL_Begin(GL_LINES);
            for(int i=node_neighbors_visible.domain.min_corner.x;i<node_neighbors_visible.domain.max_corner.x;i++) for(int j=node_neighbors_visible.domain.min_corner.y;j<node_neighbors_visible.domain.max_corner.y;j++){
                if(!node_neighbors_visible(i,j)(0)){OpenGL_Line(grid.X(TV_INT(i,j)),grid.X(TV_INT(i+1,j)));}
                if(!node_neighbors_visible(i,j)(1)){OpenGL_Line(grid.X(TV_INT(i,j)),grid.X(TV_INT(i,j+1)));}}
            OpenGL_End();

            glPopAttrib();
        }

        if(draw_density_valid_mask && density_valid_mask.Size().x) 
            opengl_density_valid_mask.Display();
        if(draw_phi_valid_mask && phi_valid_mask.Size().x) 
            opengl_phi_valid_mask.Display();
    }
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
Reinitialize(bool force)
{
    if(force || (draw && (!valid || (is_animation && (frame_loaded != frame)) || (!is_animation && (frame_loaded < 0)))))
    {
        valid = false;

        if(draw_grid_visibility){
            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(directory+"/%d/thin_shells_grid_visibility",frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,node_neighbors_visible,face_corners_visible_from_face_center_u,face_corners_visible_from_face_center_v);
        }

        if(draw_density_valid_mask){
            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(directory+"/%d/density_valid_mask",frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename)){
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,density_valid_mask);
                for(int i=density_valid_mask.domain.min_corner.x;i<density_valid_mask.domain.max_corner.x;i++) for(int j=density_valid_mask.domain.min_corner.y;j<density_valid_mask.domain.max_corner.y;j++) 
                    density_valid_mask(i,j)=!density_valid_mask(i,j); // negate
                opengl_density_valid_mask.Update();}
            else density_valid_mask.Clean_Memory();
        }

        if(draw_phi_valid_mask){
            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(directory+"/%d/phi_valid_mask",frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename)){
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,phi_valid_mask);
                for(int i=phi_valid_mask.domain.min_corner.x;i<phi_valid_mask.domain.max_corner.x;i++) for(int j=phi_valid_mask.domain.min_corner.y;j<phi_valid_mask.domain.max_corner.y;j++) 
                    phi_valid_mask(i,j)=!phi_valid_mask(i,j); // negate
                opengl_phi_valid_mask.Update();}
            else phi_valid_mask.Clean_Memory();
        }

        frame_loaded=frame;
        valid=true;
    }
}
//#####################################################################
// Function Toggle_Draw_Grid_Visibility
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
Toggle_Draw_Grid_Visibility()
{
    draw_grid_visibility=!draw_grid_visibility;
    Reinitialize(true);
}
//#####################################################################
// Function Toggle_Draw_Density_Valid_Mask
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
Toggle_Draw_Density_Valid_Mask()
{
    draw_density_valid_mask=!draw_density_valid_mask;
    Reinitialize(true);
}
//#####################################################################
// Function Toggle_Draw_Phi_Valid_Mask
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
Toggle_Draw_Phi_Valid_Mask()
{
    draw_phi_valid_mask=!draw_phi_valid_mask;
    Reinitialize(true);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T,RW>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(grid.domain);
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<float,float>;
template class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<double,double>;
}
