//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D
//##################################################################### 
#ifndef __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D__
#define __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
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

template<class T>
class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D:public OPENGL_COMPONENT<T>
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
    ARRAY<OPENGL_LEVELSET_MULTIVIEW<T>*> opengl_levelset;
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
    ARRAY<ARRAY<OPENGL_COMPONENT<T>*>,int> extra_components;
    ARRAY<VECTOR<T,3> > articulation_points;
    ARRAY<FRAME<TV> > joint_frames;
    ARRAY<PAIR<VECTOR<T,3>,VECTOR<T,3> >,int> forces_and_torques;
    VECTOR<T,3> projected_COM;
    OPENGL_SELECTION<T>* current_selection;

public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::is_animation;using OPENGL_COMPONENT<T>::stream_type;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(STREAM_TYPE stream_type,const std::string& basedir, bool use_display_lists=true);
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(STREAM_TYPE stream_type,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir, bool use_display_lists=true);
    virtual ~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D();

    bool Valid_Frame(int frame_input) const override;
    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;
    void Draw_All_Objects() override;

    virtual void Display() const override;
    bool Use_Bounding_Box() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer, int buffer_size) override;
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) override;
    void Clear_Highlight() override;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const override;

    void Turn_Smooth_Shading_On() override;
    void Turn_Smooth_Shading_Off() override;
    void Slice_Has_Changed() override;

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
    
    void Manipulate_Individual_Body_Prompt();
};

template<class T>
class OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int body_id;
    OPENGL_SELECTION<T>* body_selection;

    OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D(OPENGL_OBJECT<T>* object,const int body_id,OPENGL_SELECTION<T>* body_selection)
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::COMPONENT_RIGID_BODIES_3D,object),body_id(body_id),body_selection(body_selection)
    {}

    virtual typename OPENGL_SELECTION<T>::TYPE Actual_Type() const override 
    {return body_selection->Actual_Type();}

    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
};

template<class T>
class OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int joint_id;

    OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D(OPENGL_OBJECT<T>* object,const int joint_id)
        :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::ARTICULATED_RIGID_BODIES_JOINT_3D,object),joint_id(joint_id)
    {}

    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
};
}
#endif
