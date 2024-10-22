//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D__
#define __OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Data_Structures/PAIR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T>
class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::World_Space_Box;using OPENGL_COMPONENT<T>::viewer_callbacks;
    using OPENGL_COMPONENT<T>::viewer_dir;
    GRID<TV> grid,mac_grid,u_grid,v_grid;
    ARRAY<VECTOR<bool,2> ,VECTOR<int,2> > node_neighbors_visible;
    ARRAY<VECTOR<PAIR<bool,T>,2>,VECTOR<int,2> > face_corners_visible_from_face_center_u; // length 2, order is bottom, top
    ARRAY<VECTOR<PAIR<bool,T>,2>,VECTOR<int,2> > face_corners_visible_from_face_center_v; // length 2, order is left, right
    ARRAY<bool,VECTOR<int,2> > density_valid_mask;
    ARRAY<bool,VECTOR<int,2> > phi_valid_mask;
private:
    OPENGL_CONSTANT_COLOR_MAP<bool> invalid_color_map;
public:
    OPENGL_SCALAR_FIELD_2D<T,bool> opengl_density_valid_mask;
    OPENGL_SCALAR_FIELD_2D<T,bool> opengl_phi_valid_mask;
private:
    bool valid;
    bool draw_grid_visibility,draw_density_valid_mask,draw_phi_valid_mask;

public:
    OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D(const VIEWER_DIR& viewer_dir,GRID<TV> &grid);
    
    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;
    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid; }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    void Toggle_Draw_Grid_Visibility();
    void Toggle_Draw_Density_Valid_Mask();
    void Toggle_Draw_Phi_Valid_Mask();

private:
    void Reinitialize();
};
}
#endif
