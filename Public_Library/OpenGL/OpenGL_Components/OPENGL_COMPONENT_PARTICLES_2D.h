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
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_PARTICLES_2D(STREAM_TYPE stream_type,const std::string &filename, const std::string &filename_set_input="", bool use_ids_input = true, bool particles_stored_per_cell_uniform_input = false);
    virtual ~OPENGL_COMPONENT_PARTICLES_2D();

    bool Valid_Frame(int frame_input) const override;

    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override { return draw && valid && opengl_points->Use_Bounding_Box(); }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer, int buffer_size) override;
    virtual OPENGL_SELECTION<T>* Get_Selection_By_Id(int id,int particle_set);
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) override;
    void Clear_Highlight() override;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const override;
    OPENGL_SELECTION<T>* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION<T>* old_selection,bool& delete_selection) override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box(OPENGL_SELECTION<T>* selection) const override;
    int Get_Current_Index_Of_Selection(OPENGL_SELECTION<T>* selection) const;
    GEOMETRY_PARTICLES<TV>* Get_Particle_Set_Of_Selection(OPENGL_SELECTION<T>* selection) const;
    bool Uses_Sets() const { return use_sets; }
    void Select_Particle_By_Id(int id,int particle_set);
    void Select_Particles_By_Ids(const ARRAY<int> &ids);
    void Clear_Id_Selection();

    void Toggle_Draw_Point_Numbers();
    void Toggle_Draw_Radii();
    void Toggle_Draw_Velocities();
    void Command_Prompt();
    void Next_Set();
    void Previous_Set();
    void Toggle_Draw_Multiple_Particle_Sets();

private:
    void Reinitialize(bool force=false);
    void Apply_Id_Selection();
    ARRAY_VIEW<int>* Get_Particles_Id_Array(int set_number=0) const;

    void Command_Prompt_Response();

public:
    GEOMETRY_PARTICLES<TV>* particles;
    ARRAY<GEOMETRY_PARTICLES<TV>*> particles_multiple;
    OPENGL_POINTS_2D<T>* opengl_points;
    ARRAY<OPENGL_POINTS_2D<T>*> opengl_points_multiple;
    OPENGL_VECTOR_FIELD_2D<ARRAY<TV> > opengl_vector_field;

private:
    std::string filename;
    std::string filename_set;
    int frame_loaded;
    int set;
    int set_loaded;
    int number_of_sets;
    bool use_sets;
    bool valid;
    bool draw_velocities, have_velocities;
    bool use_ids;
    bool particles_stored_per_cell_uniform;
    bool draw_multiple_particle_sets;
    ARRAY<ARRAY<int> > selected_ids;
};

template<class T>
class OPENGL_SELECTION_COMPONENT_PARTICLES_2D:public OPENGL_SELECTION<T>
{
    typedef VECTOR<T,2> TV;
public:
    using OPENGL_SELECTION<T>::object;
    int index;  // index into particles array
    bool has_id;
    int id;
    int particle_set;
    TV location;

    OPENGL_SELECTION_COMPONENT_PARTICLES_2D(OPENGL_OBJECT<T>* object) :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::COMPONENT_PARTICLES_2D, object) {}
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

}

#endif
