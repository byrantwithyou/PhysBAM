//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_HEIGHTFIELD_1D
//#####################################################################
#ifndef __OPENGL_COMPONENT_HEIGHTFIELD_1D__
#define __OPENGL_COMPONENT_HEIGHTFIELD_1D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Grids/GRID.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T>
class OPENGL_COMPONENT_HEIGHTFIELD_1D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::World_Space_Box;
    using OPENGL_COMPONENT<T>::viewer_callbacks;using OPENGL_COMPONENT<T>::is_animation;
    OPENGL_COMPONENT_HEIGHTFIELD_1D(const GRID<TV> &grid, 
                                    const std::string& height_filename,
                                    const std::string& x_filename_input="",
                                    const std::string& ground_filename_input="",
                                    const std::string& u_filename_input="");
    virtual ~OPENGL_COMPONENT_HEIGHTFIELD_1D();

    bool Valid_Frame(int frame_input) const override;

    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;

    void Set_Scale(T scale_input);

    void Increase_Scale();
    void Decrease_Scale();
    void Increase_Displacement_Scale();
    void Decrease_Displacement_Scale();
    void Increase_Velocity_Scale();
    void Decrease_Velocity_Scale();
    void Toggle_Draw_Velocities();
    void Toggle_Draw_Points();

private:
    void Reinitialize(bool force=false);

public:
    GRID<TV> grid;
    ARRAY<T,TV_INT> *x;
    ARRAY<T,TV_INT> *ground;
    ARRAY<T,TV_INT> *u;
    ARRAY<T,TV_INT> height;

private:
    ARRAY<VECTOR<T,2> > vector_field;
    ARRAY<VECTOR<T,2> > vector_locations;
    OPENGL_VECTOR_FIELD_2D<ARRAY<VECTOR<T,2> > > opengl_vector_field;

    std::string height_filename;
    std::string x_filename;
    std::string ground_filename;
    std::string u_filename;
    T scale;
    T displacement_scale;
    int frame_loaded;
    bool valid;
    bool draw_velocities;
    bool draw_points;
    int selected_index;
};

}

#endif
