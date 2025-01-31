//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D__
#define __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D__

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
template<class TV> class DEFORMABLE_BODY_COLLECTION;

template<class T>
class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,2> TV;
protected:
    bool valid;
    int display_mode;
    bool draw_velocities;
    T velocity_scale;
    bool invalidate_deformable_objects_selection_each_frame;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_OBJECT<T>::viewer_callbacks;using OPENGL_OBJECT<T>::World_Space_Box;
    using OPENGL_COMPONENT<T>::viewer_dir;

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
    ARRAY<OPENGL_SEGMENTED_CURVE_2D<T>*> embedded_curve_objects;
    bool has_embedded_objects;
    ARRAY<ARRAY<T>*> phi_list;
    int selected_segmented_curve;
    int selected_triangulated_area;
    int selected_triangles_of_material;
    int selected_bezier_spline;
    int selected_embedded_curve;
    int selected_free_particles;
    int selected_b_spline;
    int selected_index;
    OPENGL_OBJECT<T>* selected_object;
    ARRAY<bool> active_list;

    OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D(const VIEWER_DIR& viewer_dir);
    virtual ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D();
    
    void Set_Frame() override;
    void Set_Draw(bool draw_input=true) override;

    virtual void Display() const override;
    bool Use_Bounding_Box() const override;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const override;
    virtual RANGE<VECTOR<T,3> > Selection_Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    void Clear_Selection() override;
    void Print_Selection_Info(std::ostream& output_stream) const override;
    bool Destroy_Selection_After_Frame_Change() override;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Draw_Velocities();
    void Cycle_Display_Mode();

    virtual void Reinitialize();    // Needs to be called after some state changes

private:
    void Initialize();    // Needs to be called after some state changes
};

}
#endif
