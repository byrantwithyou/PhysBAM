//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_HEIGHTFIELD_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_HEIGHTFIELD_2D__
#define __OPENGL_COMPONENT_HEIGHTFIELD_2D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T>
class OPENGL_COMPONENT_HEIGHTFIELD_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<T,2> TV2;typedef VECTOR<int,2> TV_INT2;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_HEIGHTFIELD_2D(STREAM_TYPE stream_type,const GRID<TV2> &grid, 
                                    const std::string& height_filename,
                                    const std::string& xz_filename_input="",
                                    const std::string& uv_filename_input="",
                                    int m_start_input = 0, int m_end_input = 0, int n_start_input = 0, int n_end_input = 0);
    virtual ~OPENGL_COMPONENT_HEIGHTFIELD_2D();

    bool Valid_Frame(int frame_input) const override;

    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    virtual RANGE<TV> Bounding_Box() const override;
    void Turn_Smooth_Shading_On() override;
    void Turn_Smooth_Shading_Off() override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;

    void Set_Scale(T scale_input);
    void Use_Triangle_Strip(bool use_triangle_strip_input=true);
    void Set_Grid_Filename(const std::string &grid_filename_input) {grid_filename=grid_filename_input;}

    void Increase_Scale();
    void Decrease_Scale();
    void Increase_Displacement_Scale();
    void Decrease_Displacement_Scale();
    void Increase_Velocity_Scale();
    void Decrease_Velocity_Scale();
    void Toggle_Draw_Velocities();
    void Toggle_Subdivision();

public:
    void Reinitialize(bool force=false);
private:
    void Update_Surface();

    int To_Linear_Index(int i, int j) const
    {return (i-domain.min_corner.x) + (j-domain.min_corner.y)*counts.x + 1;}

    TV_INT2 From_Linear_Index(int idx) const
    {return TV_INT2(idx % counts.x + domain.min_corner.x, (idx/counts.x)+domain.min_corner.y);}

public:
    TRIANGULATED_SURFACE<T>& triangulated_surface;
    OPENGL_TRIANGULATED_SURFACE<T> opengl_triangulated_surface;
    T vertical_offset;
    bool allow_smooth_shading;
    bool subdivide_surface;

    GRID<TV2> initial_grid;
    GRID<TV2> grid;
    ARRAY<TV2,TV_INT2 > *xz;
    ARRAY<TV2,TV_INT2 > *uv;
    ARRAY<T,TV_INT2 > height;

    ARRAY<TV> vector_field;
    ARRAY<TV> vector_locations;
    OPENGL_VECTOR_FIELD_3D<T> opengl_vector_field;

private:
    std::string grid_filename;
    std::string height_filename;
    std::string xz_filename;
    std::string uv_filename;
    T scale;
    T displacement_scale;
    int frame_loaded;
    bool valid;
    bool draw_velocities;
    RANGE<TV_INT2> domain;
    TV_INT2 counts;
    bool use_triangle_strip;
};
}

#endif
