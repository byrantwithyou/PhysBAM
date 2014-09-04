//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_PARTICLES_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_PARTICLES_3D__
#define __OPENGL_COMPONENT_PARTICLES_3D__

#include <Tools/Arrays/ARRAY.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM{

template<class T,class RW=T>
class OPENGL_COMPONENT_PARTICLES_3D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::component_name;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::World_Space_Box;
    OPENGL_COMPONENT_PARTICLES_3D(const std::string &filename, const std::string &filename_set_input="", bool use_ids_input = true, bool particles_stored_per_cell_input = false, bool particles_stored_per_cell_adaptive_input=false);
    virtual ~OPENGL_COMPONENT_PARTICLES_3D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display() const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid && opengl_points->Use_Bounding_Box(); }
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer, int buffer_size);
    virtual OPENGL_SELECTION<T>* Get_Selection_By_Id(int id);
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box(OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
    int Get_Current_Index_Of_Selection(OPENGL_SELECTION<T>* selection) const;

    void Select_Particle_By_Id(int id);
    void Select_Particles_By_Ids(const ARRAY<int> &ids);
    void Clear_Id_Selection();
    OPENGL_SELECTION<T>* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION<T>* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;

    void Toggle_Draw_Point_Numbers();
    void Toggle_Draw_Velocities();
    void Command_Prompt();
    void Set_Vector_Size(double size);
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();
    void Next_Set();
    void Previous_Set();
    void Toggle_Draw_Multiple_Particle_Sets();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D,Toggle_Draw_Point_Numbers,"Toggle draw point numbers");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D,Toggle_Draw_Velocities,"Toggle draw velocities");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D,Command_Prompt,"Command prompt");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D,Increase_Vector_Size,"Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D,Decrease_Vector_Size,"Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D,Toggle_Arrowhead,"Toggle arrow head style");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D, Next_Set, "Switch to next set");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D, Previous_Set, "Switch to previous set");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D, Toggle_Draw_Multiple_Particle_Sets, "Toggle drawing multiple particle sets");

protected:
    virtual void Reinitialize(bool force=false);
    void Apply_Id_Selection();
    ARRAY_VIEW<int>* Get_Particles_Id_Array(int set_number=0) const;

    void Command_Prompt_Response();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PARTICLES_3D, Command_Prompt_Response, "");

public:
    GEOMETRY_PARTICLES<TV>* particles;
    ARRAY<GEOMETRY_PARTICLES<TV>*> particles_multiple;
    OPENGL_POINTS_3D<T,ARRAY<VECTOR<T,3> > >* opengl_points;
    ARRAY<OPENGL_POINTS_3D<T>*> opengl_points_multiple;
    OPENGL_VECTOR_FIELD_3D<T> opengl_vector_field;

protected:
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
    bool particles_stored_per_cell_adaptive;
    bool draw_multiple_particle_sets;
    ARRAY<int> selected_ids;
};

template<class T>
class OPENGL_SELECTION_COMPONENT_PARTICLES_3D:public OPENGL_SELECTION<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_SELECTION<T>::object;
    int index;  // index into particles array
    bool has_id;
    int id;
    int particle_set;
    VECTOR<T,3> location;

    OPENGL_SELECTION_COMPONENT_PARTICLES_3D(OPENGL_OBJECT<T>* object) :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::COMPONENT_PARTICLES_3D, object) {}
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}

#endif
