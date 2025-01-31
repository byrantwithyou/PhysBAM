//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid)
    :OPENGL_COMPONENT<T>(viewer_dir,"Thin Shells Debugging"),grid(grid),
    invalid_color_map(OPENGL_COLOR::Red()),
    opengl_density_valid_mask(grid,density_valid_mask,&invalid_color_map,OPENGL_SCALAR_FIELD_3D<T,bool>::DRAW_POINTS),
    valid(false),
    draw_density_valid_mask(false),draw_node_neighbors_visible(false),draw_face_corners_visible(false)
{
    viewer_callbacks.Set("toggle_draw_grid_visibility_mode",{[this](){Toggle_Draw_Grid_Visibility_Mode();},"Toggle draw grid visibility mode"});
    viewer_callbacks.Set("toggle_draw_density_valid_mask",{[this](){Toggle_Draw_Density_Valid_Mask();},"Toggle draw density valid mask"});

    mac_grid=grid.Get_MAC_Grid();
    u_grid=grid.Get_Face_Grid(0);
    v_grid=grid.Get_Face_Grid(1);
    w_grid=grid.Get_Face_Grid(2);
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
Display() const
{
    OPENGL_COLOR node_neighbor_not_visible_color=OPENGL_COLOR::Magenta(0.5,0.5);
    OPENGL_COLOR face_corners_not_visible_from_face_center_color=OPENGL_COLOR::Magenta(1);
    if(valid && draw){
        if(draw_node_neighbors_visible || draw_face_corners_visible){
            glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_LINE_BIT | GL_CURRENT_BIT);

            if(slice && slice->Is_Slice_Mode()) slice->Enable_Clip_Planes();

            glDisable(GL_LIGHTING);

            if(draw_face_corners_visible){
                face_corners_not_visible_from_face_center_color.Send_To_GL_Pipeline();
                glLineWidth(1);
                OpenGL_Begin(GL_LINES);
                for(int i=face_corners_visible_from_face_center_u.domain.min_corner.x;i<face_corners_visible_from_face_center_u.domain.max_corner.x;i++) for(int j=face_corners_visible_from_face_center_u.domain.min_corner.y;j<face_corners_visible_from_face_center_u.domain.max_corner.y;j++) for(int k=face_corners_visible_from_face_center_u.domain.min_corner.z;k<face_corners_visible_from_face_center_u.domain.max_corner.z;k++){
                    if(!face_corners_visible_from_face_center_u(i,j,k)(0).x){OpenGL_Line(u_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j,k)));}
                    if(!face_corners_visible_from_face_center_u(i,j,k)(1).x){OpenGL_Line(u_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j+1,k)));}
                    if(!face_corners_visible_from_face_center_u(i,j,k)(2).x){OpenGL_Line(u_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j,k+1)));}
                    if(!face_corners_visible_from_face_center_u(i,j,k)(3).x){OpenGL_Line(u_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j+1,k+1)));}
                }
                for(int i=face_corners_visible_from_face_center_v.domain.min_corner.x;i<face_corners_visible_from_face_center_v.domain.max_corner.x;i++) for(int j=face_corners_visible_from_face_center_v.domain.min_corner.y;j<face_corners_visible_from_face_center_v.domain.max_corner.y;j++) for(int k=face_corners_visible_from_face_center_v.domain.min_corner.z;k<face_corners_visible_from_face_center_v.domain.max_corner.z;k++){
                    if(!face_corners_visible_from_face_center_v(i,j,k)(0).x){OpenGL_Line(v_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j,k)));}
                    if(!face_corners_visible_from_face_center_v(i,j,k)(1).x){OpenGL_Line(v_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i+1,j,k)));}
                    if(!face_corners_visible_from_face_center_v(i,j,k)(2).x){OpenGL_Line(v_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j,k+1)));}
                    if(!face_corners_visible_from_face_center_v(i,j,k)(3).x){OpenGL_Line(v_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i+1,j,k+1)));}
                }
                for(int i=face_corners_visible_from_face_center_w.domain.min_corner.x;i<face_corners_visible_from_face_center_w.domain.max_corner.x;i++) for(int j=face_corners_visible_from_face_center_w.domain.min_corner.y;j<face_corners_visible_from_face_center_w.domain.max_corner.y;j++) for(int k=face_corners_visible_from_face_center_w.domain.min_corner.z;k<face_corners_visible_from_face_center_w.domain.max_corner.z;k++){
                    if(!face_corners_visible_from_face_center_w(i,j,k)(0).x){OpenGL_Line(w_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j,k)));}
                    if(!face_corners_visible_from_face_center_w(i,j,k)(1).x){OpenGL_Line(w_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i+1,j,k)));}
                    if(!face_corners_visible_from_face_center_w(i,j,k)(2).x){OpenGL_Line(w_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j+1,k)));}
                    if(!face_corners_visible_from_face_center_w(i,j,k)(3).x){OpenGL_Line(w_grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i+1,j+1,k)));}
                }
                OpenGL_End();
            }

            if(draw_node_neighbors_visible){
                node_neighbor_not_visible_color.Send_To_GL_Pipeline();
                glLineWidth(5);
                OpenGL_Begin(GL_LINES);
                for(int i=node_neighbors_visible.domain.min_corner.x;i<node_neighbors_visible.domain.max_corner.x;i++) for(int j=node_neighbors_visible.domain.min_corner.y;j<node_neighbors_visible.domain.max_corner.y;j++) for(int k=node_neighbors_visible.domain.min_corner.z;k<node_neighbors_visible.domain.max_corner.z;k++){
                    if(!node_neighbors_visible(i,j,k)(0)){OpenGL_Line(grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i+1,j,k)));}
                    if(!node_neighbors_visible(i,j,k)(1)){OpenGL_Line(grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j+1,k)));}
                    if(!node_neighbors_visible(i,j,k)(2)){OpenGL_Line(grid.X(TV_INT(i,j,k)),grid.X(TV_INT(i,j,k+1)));}
                }
                OpenGL_End();
            }

            glPopAttrib();
        }

        if(draw_density_valid_mask && density_valid_mask.Size().x) 
            opengl_density_valid_mask.Display();
    }
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
Reinitialize()
{
    if(!draw) return;
    valid = false;
    if(draw_node_neighbors_visible || draw_face_corners_visible){
        std::string tmp_filename = viewer_dir.current_directory+"/thin_shells_grid_visibility";
        if(File_Exists(tmp_filename))
            Read_From_File(tmp_filename,node_neighbors_visible,face_corners_visible_from_face_center_u,face_corners_visible_from_face_center_v,face_corners_visible_from_face_center_w);}
    if(draw_density_valid_mask){
        std::string tmp_filename = viewer_dir.current_directory+"/density_valid_mask";
        if(File_Exists(tmp_filename)){
            Read_From_File(tmp_filename,density_valid_mask);
            for(RANGE_ITERATOR<TV::m> it(density_valid_mask.domain);it.Valid();it.Next())
                density_valid_mask(it.index)=!density_valid_mask(it.index); // negate
            opengl_density_valid_mask.Update();}
        else density_valid_mask.Clean_Memory();}
    valid=true;
}
//#####################################################################
// Function Toggle_Draw_Grid_Visibility_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
Toggle_Draw_Grid_Visibility_Mode()
{
    int mode=((int)draw_node_neighbors_visible<<1) + ((int)draw_face_corners_visible);
    mode=(mode+1)%4;
    draw_node_neighbors_visible=(mode&0x02)!=0;
    draw_face_corners_visible=(mode&0x01)!=0;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Density_Valid_Mask
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
Toggle_Draw_Density_Valid_Mask()
{
    draw_density_valid_mask=!draw_density_valid_mask;
    Reinitialize();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
Bounding_Box() const
{
    if(valid && draw) return World_Space_Box(grid.domain);
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Set_Slice
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;
    opengl_density_valid_mask.Set_Slice(slice_input);
}
//#####################################################################
// Function Slice_Has_Changed
//#####################################################################
template<class T> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>::
Slice_Has_Changed()
{
    opengl_density_valid_mask.Slice_Has_Changed();
}
namespace PhysBAM{
template class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<float>;
template class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<double>;
}
