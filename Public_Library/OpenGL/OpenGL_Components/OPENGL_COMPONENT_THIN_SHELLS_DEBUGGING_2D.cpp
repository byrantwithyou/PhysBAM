//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>::
OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D(const VIEWER_DIR& viewer_dir,GRID<TV> &grid)
    :OPENGL_COMPONENT<T>(viewer_dir,"Thin Shells Debugging"),grid(grid),
    invalid_color_map(OPENGL_COLOR::Red()),
    opengl_density_valid_mask(grid,density_valid_mask,&invalid_color_map,"density_valid",0,OPENGL_SCALAR_FIELD_2D<T,bool>::DRAW_POINTS),
    opengl_phi_valid_mask(grid,phi_valid_mask,&invalid_color_map,"phi_valid",0,OPENGL_SCALAR_FIELD_2D<T,bool>::DRAW_POINTS),
    valid(false),
    draw_grid_visibility(false),draw_density_valid_mask(false),draw_phi_valid_mask(false)
{
    viewer_callbacks.Set("toggle_draw_grid_visibility",{[this](){Toggle_Draw_Grid_Visibility();},"Toggle draw grid visibility"});
    viewer_callbacks.Set("toggle_draw_density_valid_mask",{[this](){Toggle_Draw_Density_Valid_Mask();},"Toggle draw density valid mask"});
    viewer_callbacks.Set("toggle_draw_phi_valid_mask",{[this](){Toggle_Draw_Phi_Valid_Mask();},"Toggle draw phi valid mask"});

    mac_grid=grid.Get_MAC_Grid();
    u_grid=grid.Get_Face_Grid(0);
    v_grid=grid.Get_Face_Grid(1);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>::
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
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>::
Reinitialize()
{
    if(!draw) return;
    valid = false;
    if(draw_grid_visibility){
        std::string tmp_filename = viewer_dir.current_directory+"/thin_shells_grid_visibility";
        if(File_Exists(tmp_filename))
            Read_From_File(tmp_filename,node_neighbors_visible,face_corners_visible_from_face_center_u,face_corners_visible_from_face_center_v);
    }
    if(draw_density_valid_mask){
        std::string tmp_filename = viewer_dir.current_directory+"/density_valid_mask";
        if(File_Exists(tmp_filename)){
            Read_From_File(tmp_filename,density_valid_mask);
            for(int i=density_valid_mask.domain.min_corner.x;i<density_valid_mask.domain.max_corner.x;i++) for(int j=density_valid_mask.domain.min_corner.y;j<density_valid_mask.domain.max_corner.y;j++) 
                                                                                                               density_valid_mask(i,j)=!density_valid_mask(i,j); // negate
            opengl_density_valid_mask.Update();}
        else density_valid_mask.Clean_Memory();
    }
    if(draw_phi_valid_mask){
        std::string tmp_filename = viewer_dir.current_directory+"/phi_valid_mask";
        if(File_Exists(tmp_filename)){
            Read_From_File(tmp_filename,phi_valid_mask);
            for(int i=phi_valid_mask.domain.min_corner.x;i<phi_valid_mask.domain.max_corner.x;i++) for(int j=phi_valid_mask.domain.min_corner.y;j<phi_valid_mask.domain.max_corner.y;j++) 
                                                                                                       phi_valid_mask(i,j)=!phi_valid_mask(i,j); // negate
            opengl_phi_valid_mask.Update();}
        else phi_valid_mask.Clean_Memory();
    }
    valid=true;
}
//#####################################################################
// Function Toggle_Draw_Grid_Visibility
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>::
Toggle_Draw_Grid_Visibility()
{
    draw_grid_visibility=!draw_grid_visibility;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Density_Valid_Mask
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>::
Toggle_Draw_Density_Valid_Mask()
{
    draw_density_valid_mask=!draw_density_valid_mask;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Phi_Valid_Mask
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>::
Toggle_Draw_Phi_Valid_Mask()
{
    draw_phi_valid_mask=!draw_phi_valid_mask;
    Reinitialize();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(grid.domain);
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<float>;
template class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<double>;
}
