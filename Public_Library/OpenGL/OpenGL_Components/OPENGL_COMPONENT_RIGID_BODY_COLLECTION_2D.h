//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Tamar Shinar, Rachel Weinstein, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D__
#define __OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T> class OPENGL_SEGMENTED_CURVE_2D;
template<class T> class OPENGL_TRIANGULATED_AREA;
template<class T> class OPENGL_LEVELSET_2D;
template<class T> class OPENGL_AXES;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class RIGID_BODY_COLLECTION;

template<class T>
class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D:public OPENGL_COMPONENT<T>
{
protected:
    typedef VECTOR<T,2> TV;

    std::string basedir;
    int frame_loaded;
    bool valid;
    bool show_object_names;
    bool output_positions;
    bool draw_velocity_vectors;
    bool draw_individual_axes;
    bool draw_node_velocity_vectors;
    bool draw_segmented_curve;
    bool draw_triangulated_area;
    bool draw_implicit_curve;
    bool has_init_destroy_information;
    ARRAY<int> needs_init,needs_destroy;
    bool draw_articulation_points;
    bool draw_forces_and_torques;
    bool draw_linear_muscles;
    bool need_destroy_rigid_body_collection;
    int selected_curve;
    int selected_area;
    int selected_joint_id;
    int selected_muscle_id;
    
public:
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body;
    ARRAY<int> colors;

protected:
    ARRAY<OPENGL_SEGMENTED_CURVE_2D<T>*,int> opengl_segmented_curve;
    ARRAY<OPENGL_TRIANGULATED_AREA<T>*,int> opengl_triangulated_area;
    ARRAY<OPENGL_LEVELSET_2D<T>*,int> opengl_levelset;
    ARRAY<ARRAY<OPENGL_COMPONENT<T>*>,int> extra_components;
    ARRAY<OPENGL_AXES<T>*,int> opengl_axes;
    ARRAY<bool,int> draw_object;
    ARRAY<bool,int> use_object_bounding_box;
    ARRAY<VECTOR<T,2> > positions;
    ARRAY<VECTOR<T,2> > velocity_vectors;
    ARRAY<VECTOR<T,2> > node_positions;
    ARRAY<VECTOR<T,2> > node_velocity_vectors;
    OPENGL_VECTOR_FIELD_2D<ARRAY<TV> > velocity_field;
    OPENGL_VECTOR_FIELD_2D<ARRAY<TV> > node_velocity_field;
    ARRAY<VECTOR<T,2> > articulation_points;
    ARRAY<PAIR<VECTOR<T,2>,T>,int> forces_and_torques;

public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::is_animation;
    using OPENGL_COMPONENT<T>::stream_type;using OPENGL_OBJECT<T>::viewer_callbacks;
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D(STREAM_TYPE stream_type,const std::string& basedir);
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D(STREAM_TYPE stream_type,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir);
    virtual ~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D();
    
    bool Valid_Frame(int frame_input) const override;
    void Set_Frame(int frame_input) override;
    void Set_Draw(bool draw_input = true) override;

    void Display() const override;
    bool Use_Bounding_Box() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    void Print_Selection_Info(std::ostream &output_stream) const override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;

    void Set_Draw_Object(int i, bool draw_it);  // Need to call Reinitialize after changing draw objects
    bool Get_Draw_Object(int i) const;
    void Set_Object_Color(int i, const OPENGL_COLOR &color);
    void Set_Use_Object_Bounding_Box(int i, bool use_it);
    void Set_Vector_Size(T size);

    void Toggle_Velocity_Vectors();
    void Toggle_Individual_Axes();
    void Toggle_Output_Positions();
    void Toggle_Show_Object_Names();
    void Toggle_Node_Velocity_Vectors();
    void Toggle_Draw_Mode();
    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Read_Articulated_Information(const std::string& filename);

    void Toggle_Articulation_Points();
    void Toggle_Linear_Muscles();
    void Toggle_Forces_And_Torques();

public:
    virtual void Reinitialize(const bool force=false,const bool read_geometry=true);    // Needs to be called after some state changes
protected:
    void Set_Draw_Mode(const int mode);
    int Get_Draw_Mode() const;
    void Create_Geometry(const int id);
    void Update_Geometry(const int id);
    void Destroy_Geometry(const int id);
    virtual void Update_Object_Labels();
};
}

#endif
