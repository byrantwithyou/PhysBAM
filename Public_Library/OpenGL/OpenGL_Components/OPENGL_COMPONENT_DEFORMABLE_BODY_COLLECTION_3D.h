//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D
//#####################################################################
#ifndef __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D__
#define __OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D__

#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/FREE_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_B_SPLINE_PATCH.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <OpenGL/OpenGL/OPENGL_HEXAHEDRALIZED_VOLUME.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL/OPENGL_VECTOR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D.h>
namespace PhysBAM{

template<class T> class OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D;

template<class T>
class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D:public OPENGL_COMPONENT<T>
{
    typedef VECTOR<T,3> TV;
protected:
    std::string prefix;
    int frame_loaded;
    bool valid,use_active_list,hide_unselected;
    int display_mode,display_relative_velocity_mode,number_of_segmented_curve,incremented_active_object;
    bool smooth_shading;
    int selected_vertex;
    bool invalidate_deformable_objects_selection_each_frame;
    bool own_deformable_body;
    int display_soft_bound_surface_mode,display_hard_bound_surface_mode,display_forces_mode,interaction_pair_display_mode;
public:
    using OPENGL_COMPONENT<T>::draw;using OPENGL_COMPONENT<T>::slice;using OPENGL_COMPONENT<T>::frame;
    using OPENGL_COMPONENT<T>::is_animation;using OPENGL_COMPONENT<T>::stream_type;using OPENGL_OBJECT<T>::viewer_callbacks;
    DEFORMABLE_BODY_COLLECTION<TV> &deformable_body_collection;
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>* real_selection;
    ARRAY<OPENGL_SEGMENTED_CURVE_3D<T>*> segmented_curve_objects;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> triangulated_surface_objects;
    ARRAY<OPENGL_TETRAHEDRALIZED_VOLUME<T>*> tetrahedralized_volume_objects;
    ARRAY<OPENGL_HEXAHEDRALIZED_VOLUME<T>*> hexahedralized_volume_objects;
    ARRAY<OPENGL_B_SPLINE_PATCH<T>*> b_spline_patch_objects;
    ARRAY<OPENGL_POINTS_3D<T,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >*> free_particles_objects;
    ARRAY<INDIRECT_ARRAY<ARRAY_VIEW<TV> >*> free_particles_indirect_arrays;
    bool has_tetrahedralized_volumes,has_hexahedralized_volumes;
    ARRAY<bool> active_list;
    OPENGL_COLOR_RAMP<T>* color_map_relative_velocity;
    ARRAY<TV> positions;
    ARRAY<TV> velocity_vectors;
    OPENGL_VECTOR_FIELD_3D<T> velocity_field;
    bool draw_velocity_vectors;
    OPENGL_INDEXED_COLOR_MAP *color_map;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> boundary_surface_objects;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> embedded_surface_objects;
    ARRAY<OPENGL_TRIANGULATED_SURFACE<T>*> hard_bound_boundary_surface_objects;
    bool has_embedded_objects,has_soft_bindings;
    ARRAY<REPULSION_PAIR<TV> > point_triangle_interaction_pairs;
    ARRAY<REPULSION_PAIR<TV> > edge_edge_interaction_pairs;
    ARRAY<FORCE_DATA<TV> > force_data_list;
    OPENGL_COLOR_RAMP<T>* color_map_forces;

    OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D(STREAM_TYPE stream_type,const std::string& prefix,const int start_frame);
    virtual ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D();
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Draw_All_Objects() PHYSBAM_OVERRIDE;

    void Set_Display_Modes(bool& display_triangulated_surface_objects,bool& display_tetrahedralized_volume_objects,
            bool& display_hexahedralized_volume_objects,bool& display_free_particles_objects,bool& display_boundary_surface_objects,
        bool& display_hard_bound_boundary_surface_objects) const;
    virtual void Display() const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION<T>* Get_Selection(GLuint *buffer, int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION<T>* selection) PHYSBAM_OVERRIDE;
    void Set_Selection(OPENGL_SELECTION<T>* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const PHYSBAM_OVERRIDE;
    
    void Turn_Smooth_Shading_On() PHYSBAM_OVERRIDE
    {smooth_shading=true;
    for(int i=0;i<segmented_curve_objects.m;i++)if(segmented_curve_objects(i))segmented_curve_objects(i)->Turn_Smooth_Shading_On();
    for(int i=0;i<triangulated_surface_objects.m;i++)if(triangulated_surface_objects(i))triangulated_surface_objects(i)->Turn_Smooth_Shading_On();
    for(int i=0;i<tetrahedralized_volume_objects.m;i++)if(tetrahedralized_volume_objects(i))tetrahedralized_volume_objects(i)->Turn_Smooth_Shading_On();
    for(int i=0;i<hexahedralized_volume_objects.m;i++)if(hexahedralized_volume_objects(i))hexahedralized_volume_objects(i)->Turn_Smooth_Shading_On();
    for(int i=0;i<boundary_surface_objects.m;i++)if(boundary_surface_objects(i))boundary_surface_objects(i)->Turn_Smooth_Shading_On();
    for(int i=0;i<embedded_surface_objects.m;i++)if(embedded_surface_objects(i))embedded_surface_objects(i)->Turn_Smooth_Shading_On();}
    
    void Turn_Smooth_Shading_Off() PHYSBAM_OVERRIDE
    {smooth_shading=false;
    for(int i=0;i<segmented_curve_objects.m;i++)if(segmented_curve_objects(i))segmented_curve_objects(i)->Turn_Smooth_Shading_Off();
    for(int i=0;i<triangulated_surface_objects.m;i++)if(triangulated_surface_objects(i))triangulated_surface_objects(i)->Turn_Smooth_Shading_Off();
    for(int i=0;i<tetrahedralized_volume_objects.m;i++)if(tetrahedralized_volume_objects(i))tetrahedralized_volume_objects(i)->Turn_Smooth_Shading_Off();
    for(int i=0;i<hexahedralized_volume_objects.m;i++)if(hexahedralized_volume_objects(i))hexahedralized_volume_objects(i)->Turn_Smooth_Shading_Off();
    for(int i=0;i<boundary_surface_objects.m;i++)if(boundary_surface_objects(i))boundary_surface_objects(i)->Turn_Smooth_Shading_Off();
    for(int i=0;i<embedded_surface_objects.m;i++)if(embedded_surface_objects(i))embedded_surface_objects(i)->Turn_Smooth_Shading_Off();}

    void Set_Material(const int object,const OPENGL_MATERIAL& front_material,const OPENGL_MATERIAL& back_material);
    void Set_All_Materials(const OPENGL_MATERIAL& meshfront,const OPENGL_MATERIAL& front_material,const OPENGL_MATERIAL& back_material);
    OPENGL_SELECTION<T>* Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION<T>* old_selection,bool& delete_selection) PHYSBAM_OVERRIDE;

    TRIANGULATED_SURFACE<T>& Create_Hard_Bound_Boundary_Surface(TRIANGULATED_SURFACE<T>& boundary_surface);

    void Toggle_Active_Value();
    void Toggle_Draw_Interior();
    void Toggle_Differentiate_Inverted();
    void Toggle_Draw_Subsets();
    void Toggle_Use_Active_List();
    void Toggle_Selection_Mode();
    void Toggle_Hide_Unselected();
    void Increment_Active_Object();
    void Cycle_Display_Mode();
    void Cycle_Cutaway_Mode();
    void Decrease_Cutaway_Fraction();
    void Increase_Cutaway_Fraction();
    void Create_One_Big_Triangulated_Surface_And_Write_To_File();
    void Show_Only_First();
    void Highlight_Particle();
    void Set_Vector_Size(double size);
    void Toggle_Velocity_Vectors();
    void Decrease_Vector_Size();
    void Increase_Vector_Size();
    void Update_Velocity_Field();
    void Cycle_Relative_Velocity_Mode();
    void Cycle_Forces_Mode();
    void Cycle_Hard_Bound_Surface_Display_Mode();
    void Cycle_Interaction_Pair_Display_Mode();

    virtual void Reinitialize(bool force=false,bool read_geometry=true);    // Needs to be called after some state changes
    
    void Toggle_Active_Value_Response();
    void Highlight_Particle_Response();
private:
    void Initialize();    // Needs to be called after some state changes
};

template<class T>
class OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D:public OPENGL_SELECTION<T>
{
public:
    using OPENGL_SELECTION<T>::object;
    int body_index;
    int subobject_type;
    OPENGL_OBJECT<T>* subobject;
    OPENGL_SELECTION<T>* body_selection;
    OPENGL_SELECTION<T>* saved_selection;

    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D(OPENGL_OBJECT<T>* object) :OPENGL_SELECTION<T>(OPENGL_SELECTION<T>::COMPONENT_DEFORMABLE_COLLECTION_3D,object) {}
    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    virtual typename OPENGL_SELECTION<T>::TYPE Actual_Type() const {return body_selection->Actual_Type();}
};

}
#endif
