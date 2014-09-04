//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D__
#define __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/PAIR.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_LEVELSET_MULTIVIEW.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D.h>
namespace PhysBAM{

template<class T> class OPENGL_AXES;
template<class T> class OPENGL_TRIANGULATED_SURFACE;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class RIGID_BODY_COLLECTION;

template<class T,class RW=T>
class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D:public OPENGL_COMPONENT
{
protected:
    typedef VECTOR<T,3> TV;

    std::string basedir;
    bool use_display_lists;
    int frame_loaded;
    bool valid;
    bool show_object_names;
    bool output_positions;
    bool draw_velocity_vectors;
    bool draw_angular_velocity_vectors;
    bool draw_individual_axes;
    bool draw_triangulated_surface;
    bool draw_tetrahedralized_volume;
    bool draw_implicit_surface;
    bool read_triangulated_surface,read_implicit_surface,read_tetrahedralized_volume;
    ARRAY<int> needs_init,needs_destroy;
    bool has_init_destroy_information;

public:
    RIGID_BODY_COLLECTION<TV> &rigid_body_collection;
    ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body;
    ARRAY<OPENGL_COLOR,int> opengl_colors;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> opengl_triangulated_surface;
protected:
    ARRAY<OPENGL_TETRAHEDRALIZED_VOLUME<T>*> opengl_tetrahedralized_volume;
    ARRAY<OPENGL_LEVELSET_MULTIVIEW<T,RW>*> opengl_levelset;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> opengl_octree_levelset_surface;
    ARRAY<OPENGL_AXES<T>*> opengl_axes;
    ARRAY<bool> draw_object;
    ARRAY<bool> use_object_bounding_box;
    ARRAY<VECTOR<T,3> > positions;
    ARRAY<VECTOR<T,3> > velocity_vectors;
    ARRAY<VECTOR<T,3> > angular_velocity_vectors;
    OPENGL_VECTOR_FIELD_3D<T> velocity_field;
    OPENGL_VECTOR_FIELD_3D<T> angular_velocity_field;
    bool need_destroy_rigid_body_collection;
    bool one_sided;
    OPENGL_INDEXED_COLOR_MAP *front_color_map,*back_color_map;
    bool draw_simplicial_object_particles;
    bool draw_articulation_points;
    int draw_joint_frames;
    bool draw_forces_and_torques;
    ARRAY<ARRAY<OPENGL_COMPONENT*>,int> extra_components;
    ARRAY<VECTOR<T,3> > articulation_points;
    ARRAY<FRAME<TV> > joint_frames;
    ARRAY<PAIR<VECTOR<T,3>,VECTOR<T,3> >,int> forces_and_torques;
    VECTOR<T,3> projected_COM;
    OPENGL_SELECTION* current_selection;

public:
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(const std::string& basedir, bool use_display_lists=true);
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir, bool use_display_lists=true);
    virtual ~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Draw_All_Objects() PHYSBAM_OVERRIDE;

    virtual void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION* Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION* selection) const PHYSBAM_OVERRIDE;

    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE;
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE;
    void Slice_Has_Changed() PHYSBAM_OVERRIDE;

    void Read_Hints(const std::string& filename);

    void Set_Draw_Object(int i, bool draw_it);  // Need to call Reinitialize after changing draw objects
    bool Get_Draw_Object(int i) const;
    void Set_Object_Material(int i, const OPENGL_MATERIAL &front_material_input);
    void Set_Object_Material(int i, const OPENGL_MATERIAL &front_material_input, const OPENGL_MATERIAL &back_material_input);
    void Set_Use_Object_Bounding_Box(int i, bool use_it);
    void Set_Vector_Size(double size);

    void Toggle_Velocity_Vectors();
    void Toggle_Angular_Velocity_Vectors();
    void Toggle_Individual_Axes();
    void Toggle_Output_Positions();
    void Toggle_Show_Object_Names();
    void Turn_Off_Individual_Smooth_Shading();
    void Manipulate_Individual_Body();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Draw_Values();
    void Toggle_One_Sided();
    void Toggle_Draw_Particles();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Velocity_Vectors, "Toggle velocity vectors");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Angular_Velocity_Vectors, "Toggle angular velocity vectors");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Individual_Axes, "Toggle individual axes");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Output_Positions, "Toggle output positions");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Show_Object_Names, "Toggle show object names");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Turn_Off_Individual_Smooth_Shading, "Turn off individual smooth shading");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Manipulate_Individual_Body, "Manipulate individual body");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Draw_Mode, "Toggle draw mode");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Draw_Values, "Toggle draw values");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_One_Sided, "Toggle one/two sided drawing");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Draw_Particles, "Toggle drawing of simplicial object's vertices");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Articulation_Points, "Toggle articulation points");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Joint_Frames, "Toggle joint frames");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Toggle_Forces_And_Torques, "Toggle forces and torques");

    void Read_Articulated_Information(const std::string& filename);
    void Update_Articulation_Points();

    void Toggle_Articulation_Points();
    void Toggle_Joint_Frames();
    void Toggle_Forces_And_Torques();

    void Resize_Structures(const int size);
    virtual void Reinitialize(const bool force=false,const bool read_geometry=true);    // Needs to be called after some state changes
    virtual void Reinitialize_Without_Files(const bool force=false);    // Needs to be called after some state changes
    void Update_Bodies(const bool update_arb_points=true);
    void Initialize_One_Body(const int body_id,const bool force=false);
protected:
    void Initialize();
    void Set_Draw_Mode(const int mode);
    int Get_Draw_Mode() const;
    void Create_Geometry(const int id);
    void Update_Geometry(const int id);
    void Destroy_Geometry(const int id);
    void Update_Object_Labels();
    void Initialize_Display_Lists();
    void Toggle_Draw_Mode();

    void Turn_Off_Individual_Smooth_Shading_Prompt();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Turn_Off_Individual_Smooth_Shading_Prompt, "");
    
    void Manipulate_Individual_Body_Prompt();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D, Manipulate_Individual_Body_Prompt, "");
};

template<class T>
class OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D:public OPENGL_SELECTION
{
public:
    int body_id;
    OPENGL_SELECTION* body_selection;

    OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D(OPENGL_OBJECT* object,const int body_id,OPENGL_SELECTION* body_selection)
        :OPENGL_SELECTION(OPENGL_SELECTION::COMPONENT_RIGID_BODIES_3D,object),body_id(body_id),body_selection(body_selection)
    {}

    virtual OPENGL_SELECTION::TYPE Actual_Type() const 
    {return body_selection->Actual_Type();}

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

template<class T>
class OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D:public OPENGL_SELECTION
{
public:
    int joint_id;

    OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D(OPENGL_OBJECT* object,const int joint_id)
        :OPENGL_SELECTION(OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_3D,object),joint_id(joint_id)
    {}

    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};
}
#endif
