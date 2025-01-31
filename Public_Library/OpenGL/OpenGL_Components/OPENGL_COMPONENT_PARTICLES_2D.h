//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_PARTICLES_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_PARTICLES_2D__
#define __OPENGL_COMPONENT_PARTICLES_2D__

#include <Core/Arrays/ARRAY.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_POINTS_2D.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{

template<class T>
class OPENGL_COMPONENT_PARTICLES_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::World_Space_Box;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::viewer_dir;
    using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_PARTICLES_2D(const VIEWER_DIR& viewer_dir, const std::string &filename_input="", const std::string &filename_set_input="", bool use_ids_input = true, bool particles_stored_per_cell_uniform_input = false);
    virtual ~OPENGL_COMPONENT_PARTICLES_2D();

    void Set_Frame() override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid && opengl_points->Use_Bounding_Box(); }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;
    void Clear_Selection() override;
    void Print_Selection_Info(std::ostream &output_stream) const override;
    bool Destroy_Selection_After_Frame_Change() override;

    void Toggle_Draw_Point_Numbers();
    void Toggle_Draw_Radii();
    void Toggle_Draw_Velocities();
    void Command_Prompt();
    void Next_Set();
    void Previous_Set();
    void Toggle_Draw_Multiple_Particle_Sets();

protected:
    virtual void Reinitialize();
    ARRAY_VIEW<int>* Get_Particles_Id_Array(int set_number=0) const;

    void Command_Prompt_Response();

public:
    GEOMETRY_PARTICLES<TV>* particles;
    ARRAY<GEOMETRY_PARTICLES<TV>*> particles_multiple;
    OPENGL_POINTS_2D<T>* opengl_points;
    ARRAY<OPENGL_POINTS_2D<T>*> opengl_points_multiple;
    OPENGL_VECTOR_FIELD_2D<ARRAY<TV> > opengl_vector_field;

    std::string filename;
    std::string filename_set;
    int set;
    int set_loaded;
    int number_of_sets;
    bool use_sets;
    bool valid;
    bool draw_velocities, have_velocities;
    bool use_ids;
    bool particles_stored_per_cell_uniform;
    bool draw_multiple_particle_sets;
    int selected_set;
};
}

#endif
