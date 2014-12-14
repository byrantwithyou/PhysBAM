//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D__
#define __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D__

#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Deformables/Particles/FREE_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_POINTS_2D.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T> class OPENGL_SEGMENTED_CURVE_2D;
template<class T,int d> class OPENGL_BEZIER_SPLINE_2D;
template<class T,int d> class OPENGL_B_SPLINE_2D;
template<class T> class OPENGL_TRIANGULATED_AREA;
template<class TV> class OPENGL_ADAPTIVE_NODE_SCALAR_FIELD;
template<class TV> class DEFORMABLE_BODY_COLLECTION;

template<class T,class RW=T>
class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;
protected:
    std::string prefix;
    int frame_loaded;
    bool valid;
    int display_mode;
    bool draw_velocities;
    T velocity_scale;
    bool invalidate_deformable_objects_selection_each_frame;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;using OPENGL_COMPONENT<T>::is_animation;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    ARRAY<OPENGL_SEGMENTED_CURVE_2D<T>*> segmented_curve_objects;
    ARRAY<OPENGL_BEZIER_SPLINE_2D<T,3>*> bezier_spline_objects;
    ARRAY<OPENGL_B_SPLINE_2D<T,3>*> b_spline_objects;
    ARRAY<OPENGL_TRIANGULATED_AREA<T>*> triangulated_area_objects;
    ARRAY<OPENGL_TRIANGULATED_AREA<T>*> triangles_of_material_objects;
    ARRAY<OPENGL_POINTS_2D<T,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >*> free_particles_objects;
    ARRAY<INDIRECT_ARRAY<ARRAY_VIEW<TV> >*> free_particles_indirect_arrays;
    OPENGL_VECTOR_FIELD_2D<ARRAY_VIEW<TV> > velocity_field;
    OPENGL_INDEXED_COLOR_MAP *color_map;
    COLLISION_BODY_COLLECTION<TV>& collision_body_list;
    ARRAY<OPENGL_SEGMENTED_CURVE_2D<T>*> embedded_curve_objects;
    bool has_embedded_objects;
    ARRAY<ARRAY<T>*> phi_list;

    OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D(const std::string& prefix,const int start_frame);
    virtual ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D();
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;

    virtual void Display() const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint* buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
    OPENGL_SELECTION<T>* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION<T>* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;

    void Set_Vector_Size(const T vector_size);

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Draw_Velocities();
    void Cycle_Display_Mode();
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D,Increase_Vector_Size,"Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D,Decrease_Vector_Size,"Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D,Toggle_Draw_Velocities,"Toggle draw velocities");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D,Cycle_Display_Mode,"Cycle embedded display mode");

    virtual void Reinitialize(bool force=false);    // Needs to be called after some state changes

private:
    void Initialize();    // Needs to be called after some state changes
};

template<class T>
class OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_2D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int body_index;
    int subobject_type;
    OPENGL_OBJECT<T>* subobject;
    OPENGL_SELECTION<T>* body_selection;

    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_2D(OPENGL_OBJECT<T>* object):OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::COMPONENT_DEFORMABLE_OBJECT_2D, object){}
    RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    typename OPENGL_SELECTION<T>::TYPE Actual_Type() const PHYSBAM_OVERRIDE {return body_selection->Actual_Type();}
};

}
#endif
