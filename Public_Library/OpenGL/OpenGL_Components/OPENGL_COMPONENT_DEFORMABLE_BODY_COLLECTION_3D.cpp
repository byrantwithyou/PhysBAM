//#####################################################################
// Copyright 2004-2009, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_FREE_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D
//#####################################################################
template<class T> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D(STREAM_TYPE stream_type,const std::string& prefix,const int start_frame)
    :OPENGL_COMPONENT<T>(stream_type,"Deformable Object List"),prefix(prefix),frame_loaded(-1),valid(false),use_active_list(false),hide_unselected(false),display_mode(0),display_relative_velocity_mode(0),number_of_segmented_curve(1),
    incremented_active_object(0),smooth_shading(false),selected_vertex(-1),
    display_hard_bound_surface_mode(0),display_forces_mode(0),interaction_pair_display_mode(0),
    collision_body_list(*new COLLISION_BODY_COLLECTION<TV>),
    deformable_body_collection(*new DEFORMABLE_BODY_COLLECTION<TV>(collision_body_list)),real_selection(0),
    has_tetrahedralized_volumes(false),has_hexahedralized_volumes(false),
    velocity_field(stream_type,velocity_vectors,positions,OPENGL_COLOR::Cyan(),.25,false,false),
    color_map(OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map()),
    has_embedded_objects(false),has_soft_bindings(false)
{
    viewer_callbacks.Set("toggle_active_value",{[this](){Toggle_Active_Value();},"Toggle viewing of elements"});
    viewer_callbacks.Set("toggle_draw_interior",{[this](){Toggle_Draw_Interior();},"Toggle view of interior elements for tetraheralized volumes"});
    viewer_callbacks.Set("toggle_differentiate_inverted",{[this](){Toggle_Differentiate_Inverted();},"Toggle use of different color for inverted tetrahera"});
    viewer_callbacks.Set("toggle_draw_subsets",{[this](){Toggle_Draw_Subsets();},"Toggle drawing of subset tets and particles for tetraheralized volumes"});
    viewer_callbacks.Set("toggle_hide_unselected",{[this](){Toggle_Hide_Unselected();},"Toggle drawing of the selected regions"});
    viewer_callbacks.Set("toggle_use_active_list",{[this](){Toggle_Use_Active_List();},"Toggle drawing subset of the deformable objects in the list"});
    viewer_callbacks.Set("toggle_selection_mode",{[this](){Toggle_Selection_Mode();},"Toggle selecting a whole segment or just one part"});
    viewer_callbacks.Set("increment_active_object",{[this](){Increment_Active_Object();},"Increment deformable object being drawn"});
    viewer_callbacks.Set("cycle_display_mode",{[this](){Cycle_Display_Mode();},"Cycle embedded display mode"});
    viewer_callbacks.Set("cycle_cutaway_mode",{[this](){Cycle_Cutaway_Mode();},"Cycle cutaway mode"});
    viewer_callbacks.Set("decrease_cutaway_fraction",{[this](){Decrease_Cutaway_Fraction();},"Decrease_Cutaway_Fraction"});
    viewer_callbacks.Set("increase_cutaway_fraction",{[this](){Increase_Cutaway_Fraction();},"Increase_Cutaway_Fraction"});
    viewer_callbacks.Set("create_one_big_triangulated_surface_and_write_to_file",{[this](){Create_One_Big_Triangulated_Surface_And_Write_To_File();},"Make one big boundary tri surface"});
    viewer_callbacks.Set("cycle_relative_velocity_mode",{[this](){Cycle_Relative_Velocity_Mode();},"Visualize relative velocity"});
    viewer_callbacks.Set("show_only_first",{[this](){Show_Only_First();},"Show only first deformable object"});
    viewer_callbacks.Set("highlight_particle",{[this](){Highlight_Particle();},"Highlight a particle"});
    viewer_callbacks.Set("toggle_velocity_vectors",{[this](){Toggle_Velocity_Vectors();},"Toggle particle velocity vectors"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("cycle_forces_mode",{[this](){Cycle_Forces_Mode();},"Cycle forces display mode"});
    viewer_callbacks.Set("cycle_hard_bound_surface_display_mode",{[this](){Cycle_Hard_Bound_Surface_Display_Mode();},"Cycle hard bound embedded surface display mode"});
    viewer_callbacks.Set("cycle_interaction_pair_display_mode",{[this](){Cycle_Interaction_Pair_Display_Mode();},"Cycle display of interaction pairs"});

    // check for per frame particles
    if(FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/deformable_object_structures",prefix.c_str(),start_frame)))
        invalidate_deformable_objects_selection_each_frame=true;
    else invalidate_deformable_objects_selection_each_frame=false;

    color_map_relative_velocity=new OPENGL_COLOR_RAMP<T>();
    color_map_relative_velocity->Add_Color(-.05,OPENGL_COLOR(1,0,0));
    color_map_relative_velocity->Add_Color(0,OPENGL_COLOR(0.5,0.5,0.5));
    color_map_relative_velocity->Add_Color(.05,OPENGL_COLOR(0,0,1));

    is_animation=true;
    color_map_forces=new OPENGL_COLOR_RAMP<T>();
    color_map_forces->Add_Color(-1.,OPENGL_COLOR(0,0,1));
    color_map_forces->Add_Color(-.01,OPENGL_COLOR(0,1,1));
    color_map_forces->Add_Color(.0,OPENGL_COLOR(0,0,0));
    color_map_forces->Add_Color(.01,OPENGL_COLOR(1,1,0));
    color_map_forces->Add_Color(1.,OPENGL_COLOR(1,0,0));
    Initialize();
}
//#####################################################################
// Function ~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D
//#####################################################################
template<class T> OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
~OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D()
{
    delete color_map;
    delete color_map_relative_velocity;
    segmented_curve_objects.Delete_Pointers_And_Clean_Memory();
    triangulated_surface_objects.Delete_Pointers_And_Clean_Memory();
    tetrahedralized_volume_objects.Delete_Pointers_And_Clean_Memory();
    hexahedralized_volume_objects.Delete_Pointers_And_Clean_Memory();
    free_particles_objects.Delete_Pointers_And_Clean_Memory();
    free_particles_indirect_arrays.Delete_Pointers_And_Clean_Memory();
    delete color_map_forces;
    delete &deformable_body_collection;
    delete &collision_body_list;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Initialize()
{
    draw_velocity_vectors=false;
}
//#####################################################################
// Function Set_Material
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Set_Material(const int object,const OPENGL_MATERIAL& front_material,const OPENGL_MATERIAL& back_material)
{
    if(triangulated_surface_objects(object)){
        triangulated_surface_objects(object)->front_material=front_material;
        triangulated_surface_objects(object)->back_material=back_material;
        triangulated_surface_objects(object)->Set_Two_Sided(true);}
    else if(tetrahedralized_volume_objects(object))
        tetrahedralized_volume_objects(object)->material=front_material;
    else if(hexahedralized_volume_objects(object))
        hexahedralized_volume_objects(object)->material=front_material;
}
//#####################################################################
// Function Set_All_Materials
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Set_All_Materials(const OPENGL_MATERIAL& meshfront,const OPENGL_MATERIAL& front_material,const OPENGL_MATERIAL& back_material)
{
    for(int i=0;i<deformable_body_collection.structures.m;i++) Set_Material(i,front_material,back_material);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Reinitialize(bool force,bool read_geometry)
{
    if(!(draw && (force || (is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded<0)))) return;
    static bool first_time=true;
    std::string frame_string=LOG::sprintf("%s/%d/",prefix.c_str(),frame);
    std::string static_frame_string=frame_string;
    int static_frame=FILE_UTILITIES::File_Exists(frame_string+"deformable_object_structures")?frame:-1;
    bool read_static_variables=static_frame!=-1 || first_time || !deformable_body_collection.structures.m;
    if(read_geometry) deformable_body_collection.Read(stream_type,prefix,prefix,frame,static_frame,read_static_variables,true);
    if(FILE_UTILITIES::File_Exists(frame_string+"/interaction_pairs") && interaction_pair_display_mode)
        FILE_UTILITIES::Read_From_File(stream_type,frame_string+"/interaction_pairs",point_triangle_interaction_pairs,edge_edge_interaction_pairs);
    else{
        point_triangle_interaction_pairs.Remove_All();edge_edge_interaction_pairs.Remove_All();}
    
    std::string filename=frame_string+"/deformable_object_force_data";
    if(FILE_UTILITIES::File_Exists(filename)){
        if(first_time) LOG::cout<<"reading "<<filename<<std::endl;
        FILE_UTILITIES::Read_From_File(stream_type,filename,force_data_list);}
    else force_data_list.Remove_All();

    if(read_static_variables){
        int m=deformable_body_collection.structures.m;active_list.Resize(m,true,true,true);
        segmented_curve_objects.Delete_Pointers_And_Clean_Memory();segmented_curve_objects.Resize(m);
        triangulated_surface_objects.Delete_Pointers_And_Clean_Memory();triangulated_surface_objects.Resize(m);
        tetrahedralized_volume_objects.Delete_Pointers_And_Clean_Memory();tetrahedralized_volume_objects.Resize(m);
        hexahedralized_volume_objects.Delete_Pointers_And_Clean_Memory();hexahedralized_volume_objects.Resize(m);
        free_particles_objects.Delete_Pointers_And_Clean_Memory();free_particles_objects.Resize(m);
        free_particles_indirect_arrays.Delete_Pointers_And_Clean_Memory();free_particles_indirect_arrays.Resize(m);
        embedded_surface_objects.Delete_Pointers_And_Clean_Memory();embedded_surface_objects.Resize(m);
        boundary_surface_objects.Delete_Pointers_And_Clean_Memory();boundary_surface_objects.Resize(m);
        hard_bound_boundary_surface_objects.Delete_Pointers_And_Clean_Memory();hard_bound_boundary_surface_objects.Resize(m);
        int color_map_index=15;
        for(int i=0;i<m;i++){
            STRUCTURE<TV>* structure=deformable_body_collection.structures(i);
            if(SEGMENTED_CURVE<TV>* segmented_curve=dynamic_cast<SEGMENTED_CURVE<TV>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": segmented curve\n";
                segmented_curve_objects(i)=new OPENGL_SEGMENTED_CURVE_3D<T>(stream_type,*segmented_curve,OPENGL_COLOR((T).5,(T).25,0));
                segmented_curve_objects(i)->use_solid_color=false;}
            else if(TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(structure)){
                if(first_time && triangulated_surface->mesh.elements.m)
                    LOG::cout<<"object "<<i<<": triangulated surface, range = "<<triangulated_surface->mesh.elements.Flattened().Min()<<" "<<triangulated_surface->mesh.elements.Flattened().Max()<<"\n";
                else LOG::cout<<"object "<<i<<": triangulated surface, empty\n";
                triangulated_surface->mesh.Initialize_Segment_Mesh();
                ARRAY<OPENGL_COLOR> front_colors,back_colors;
                front_colors.Append(OPENGL_COLOR::Yellow());back_colors.Append(OPENGL_COLOR::Magenta());
                front_colors.Append(OPENGL_COLOR::Green());back_colors.Append(OPENGL_COLOR::Cyan());
                triangulated_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,*triangulated_surface,false,
                    OPENGL_MATERIAL::Metal(front_colors(i%front_colors.Size())),OPENGL_MATERIAL::Metal(back_colors(i%back_colors.Size())));}
            else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(structure)){
                if(first_time && tetrahedralized_volume->mesh.elements.m)
                    LOG::cout<<"object "<<i<<": tetrahedralized_volume, range = "<<tetrahedralized_volume->mesh.elements.Flattened().Min()<<" "<<tetrahedralized_volume->mesh.elements.Flattened().Max()<<"\n";
                else LOG::cout<<"object "<<i<<": tetrahedralized_volume, empty\n";
                tetrahedralized_volume_objects(i)=new OPENGL_TETRAHEDRALIZED_VOLUME<T>(stream_type,&tetrahedralized_volume->mesh,&(deformable_body_collection.particles),
                    OPENGL_MATERIAL::Metal(OPENGL_COLOR::Red(.7f)),OPENGL_MATERIAL::Metal(OPENGL_COLOR::Green(.7f)));}
            else if(HEXAHEDRALIZED_VOLUME<T>* hexahedralized_volume=dynamic_cast<HEXAHEDRALIZED_VOLUME<T>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": hexahedralized_volume\n";
                hexahedralized_volume_objects(i)=new OPENGL_HEXAHEDRALIZED_VOLUME<T>(stream_type,&hexahedralized_volume->mesh,&(deformable_body_collection.particles),
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()),OPENGL_MATERIAL::Matte(OPENGL_COLOR::Green()));}
            else if(FREE_PARTICLES<TV>* fp=dynamic_cast<FREE_PARTICLES<TV>*>(structure)){
                free_particles_indirect_arrays(i)=new INDIRECT_ARRAY<ARRAY_VIEW<TV> >(deformable_body_collection.particles.X,fp->nodes);
                free_particles_objects(i)=new OPENGL_FREE_PARTICLES<TV>(stream_type,deformable_body_collection,*free_particles_indirect_arrays(i),color_map->Lookup(color_map_index--));}
            if(EMBEDDED_MATERIAL_SURFACE<TV,2>* embedding=dynamic_cast<EMBEDDED_MATERIAL_SURFACE<TV,2>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": embedded triangulated surface\n";
                boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,embedding->material_surface,false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Yellow()),OPENGL_MATERIAL::Matte(OPENGL_COLOR::Cyan()));
                hard_bound_boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,Create_Hard_Bound_Boundary_Surface(embedding->material_surface),false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Magenta(.5f)));
                embedding->embedded_object.simplicial_object.mesh.Initialize_Segment_Mesh();
                triangulated_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,embedding->embedded_object.simplicial_object,false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()),OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()));}
            else if(EMBEDDED_MATERIAL_SURFACE<TV,3>* embedding=dynamic_cast<EMBEDDED_MATERIAL_SURFACE<TV,3>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": embedded tetrahedralized volume\n";
                boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,embedding->material_surface,false,OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Cyan()));
                embedded_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,embedding->embedded_object.embedded_object,false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()),OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()));
                hard_bound_boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,Create_Hard_Bound_Boundary_Surface(embedding->material_surface),false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Magenta(.5f)));
                embedding->embedded_object.simplicial_object.mesh.Initialize_Neighbor_Nodes();
                tetrahedralized_volume_objects(i)=new OPENGL_TETRAHEDRALIZED_VOLUME<T>(stream_type,&embedding->embedded_object.simplicial_object.mesh,
                    &deformable_body_collection.particles,OPENGL_MATERIAL::Matte(OPENGL_COLOR::Red()),OPENGL_MATERIAL::Matte(OPENGL_COLOR::Green()));}
            else if(EMBEDDING<TV>* embedding=dynamic_cast<EMBEDDING<TV>*>(structure)){
                if(first_time) LOG::cout<<"object "<<i<<": embedding\n";
                boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,embedding->material_surface,false,OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Cyan()));
                hard_bound_boundary_surface_objects(i)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,Create_Hard_Bound_Boundary_Surface(embedding->material_surface),false,
                    OPENGL_MATERIAL::Matte(OPENGL_COLOR::Magenta(.5f)));}
            else{if(first_time) LOG::cout<<"object "<<i<<": object unrecognized\n";}}}

    for(int i=0;i<deformable_body_collection.structures.m;i++){
        std::string i_string=LOG::sprintf("%d",i);
        if(tetrahedralized_volume_objects(i)) has_tetrahedralized_volumes=true;
        if(hexahedralized_volume_objects(i)) has_hexahedralized_volumes=true;
        if(boundary_surface_objects(i)) has_embedded_objects=true;
#if 0 // TODO: Fix me
        if(deformable_body_collection.soft_bindings.bindings.m) has_soft_bindings=true;
#endif
        if(tetrahedralized_volume_objects(i)){
            std::string filename=LOG::sprintf("%s/%d/subset_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File(stream_type,filename,tetrahedralized_volume_objects(i)->subset);
            filename=LOG::sprintf("%s/%d/colliding_nodes_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File(stream_type,filename,tetrahedralized_volume_objects(i)->subset_particles);}
        else if(hexahedralized_volume_objects(i)){
            std::string filename=LOG::sprintf("%s/%d/subset_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File(stream_type,filename,hexahedralized_volume_objects(i)->subset);
            filename=LOG::sprintf("%s/%d/colliding_nodes_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File(stream_type,filename,hexahedralized_volume_objects(i)->subset_particles);
            filename=LOG::sprintf("%s/%d/directions_%d",prefix.c_str(),frame,i);
            if(FILE_UTILITIES::File_Exists(filename))FILE_UTILITIES::Read_From_File(stream_type,filename,hexahedralized_volume_objects(i)->vectors_at_hex_centers);}}
    if(smooth_shading){
        for(int i=0;i<triangulated_surface_objects.m;i++) if(triangulated_surface_objects(i)) triangulated_surface_objects(i)->Initialize_Vertex_Normals();
        for(int i=0;i<tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i))tetrahedralized_volume_objects(i)->Initialize_Vertex_Normals();
        for(int i=0;i<boundary_surface_objects.m;i++) if(boundary_surface_objects(i))boundary_surface_objects(i)->Initialize_Vertex_Normals();
        for(int i=0;i<embedded_surface_objects.m;i++) if(embedded_surface_objects(i))embedded_surface_objects(i)->Initialize_Vertex_Normals();}
    for(int i=0;i<free_particles_indirect_arrays.m;i++) if(free_particles_indirect_arrays(i)){
        ARRAY_VIEW<TV> tmp(deformable_body_collection.particles.X);
        free_particles_indirect_arrays(i)->array.Exchange(tmp);}

    Update_Velocity_Field();
    frame_loaded=frame;
    valid=true;
    first_time=false;

    for(number_of_segmented_curve=0;deformable_body_collection.template Find_Structure<SEGMENTED_CURVE<TV>*>(number_of_segmented_curve);number_of_segmented_curve++);
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/deformable_object_particles",prefix.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    if(draw_input) active_list.Fill(true);
    Update_Velocity_Field();
    Reinitialize();
}
//#####################################################################
// Function Draw_All_Objects
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Draw_All_Objects()
{
    Set_Draw(true);
}
//#####################################################################
// Function Set_Display_Modes
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Set_Display_Modes(bool& display_triangulated_surface_objects,bool& display_tetrahedralized_volume_objects,
        bool& display_hexahedralized_volume_objects,bool& display_free_particles_objects,bool& display_boundary_surface_objects,
        bool& display_hard_bound_boundary_surface_objects) const
{
    display_triangulated_surface_objects=display_mode!=1;
    display_tetrahedralized_volume_objects=display_mode!=2;
    display_hexahedralized_volume_objects=display_mode!=2;
    display_free_particles_objects=display_mode!=1;
    if(has_embedded_objects || has_soft_bindings){
        display_triangulated_surface_objects=display_mode!=1 && display_soft_bound_surface_mode!=2;
        display_tetrahedralized_volume_objects=display_mode!=1 && display_mode!=3;// && display_soft_bound_surface_mode;
        display_hexahedralized_volume_objects=display_mode!=1;
        display_free_particles_objects=display_mode!=1;}
    else{
        display_triangulated_surface_objects=display_mode!=1;
        display_tetrahedralized_volume_objects=display_mode!=2;
        display_hexahedralized_volume_objects=display_mode!=2;
        display_free_particles_objects=display_mode!=1;}
    if(has_embedded_objects || has_soft_bindings){
        display_boundary_surface_objects=display_mode!=2 && display_hard_bound_surface_mode!=2;
        display_hard_bound_boundary_surface_objects=display_mode!=2 && display_hard_bound_surface_mode;}
    else{
        display_boundary_surface_objects=display_mode!=1;
        display_hard_bound_boundary_surface_objects=display_mode!=1 && display_hard_bound_surface_mode;}
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Display() const
{
    if(!draw || !valid) return;
    if(slice && slice->Is_Slice_Mode()){
        glPushAttrib(GL_ENABLE_BIT);
        slice->Enable_Clip_Planes();}

    bool display_triangulated_surface_objects,display_tetrahedralized_volume_objects,display_hexahedralized_volume_objects,
         display_boundary_surface_objects,display_hard_bound_boundary_surface_objects,display_free_particles_objects;
    Set_Display_Modes(display_triangulated_surface_objects,display_tetrahedralized_volume_objects,
            display_hexahedralized_volume_objects,display_boundary_surface_objects,display_hard_bound_boundary_surface_objects,display_free_particles_objects);

    for(int i=0;i<segmented_curve_objects.m;i++){
        if(!active_list(i)) continue;
        glPushName(i);
        if(segmented_curve_objects(i)){
            glPushName(1);
            segmented_curve_objects(i)->hide_unselected=hide_unselected;
            if(real_selection) segmented_curve_objects(i)->parent_curve=dynamic_cast<OPENGL_SEGMENTED_CURVE_3D<T>*>(real_selection->subobject);
            else segmented_curve_objects(i)->parent_curve=0;
            segmented_curve_objects(i)->Display();
            glPopName();}
        if(triangulated_surface_objects(i) && display_triangulated_surface_objects){
            glPushName(2);triangulated_surface_objects(i)->Display();glPopName();}
        if(tetrahedralized_volume_objects(i) && display_tetrahedralized_volume_objects){
            glPushName(3);tetrahedralized_volume_objects(i)->Display();glPopName();}
        if(hexahedralized_volume_objects(i) && display_hexahedralized_volume_objects){
            glPushName(3);hexahedralized_volume_objects(i)->Display();glPopName();}
        if(free_particles_objects(i) && display_free_particles_objects){glPushName(6);free_particles_objects(i)->Display();glPopName();}
        glPopName();}

    if(slice && slice->Is_Slice_Mode()) glPopAttrib();

    if(selected_vertex>=0) OPENGL_SELECTION<T>::Draw_Highlighted_Vertex(deformable_body_collection.particles.X(selected_vertex),selected_vertex);

    if(draw_velocity_vectors) velocity_field.Display();

    // Visualize relative velocity on edges
    if(display_relative_velocity_mode){
        glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
        glEnable(GL_BLEND);
        glDisable(GL_LIGHTING);
        OpenGL_Begin(GL_LINES);
        SEGMENTED_CURVE<TV>* segmented_curve=deformable_body_collection.template Find_Structure<SEGMENTED_CURVE<TV>*>(display_relative_velocity_mode);
        for(int j=0;j<segmented_curve->mesh.elements.m;j++){int p1=segmented_curve->mesh.elements(j)(0),p2=segmented_curve->mesh.elements(j)(1);
            TV relative_velocity=deformable_body_collection.particles.V(p2)-deformable_body_collection.particles.V(p1);
            TV edge_vector=deformable_body_collection.particles.X(p2)-deformable_body_collection.particles.X(p1);
            T edge_length=edge_vector.Magnitude();
            if(edge_length>(T)0){
                OPENGL_COLOR edge_color=color_map_relative_velocity->Lookup(TV::Dot_Product(relative_velocity,edge_vector)/edge_length);
                edge_color.Send_To_GL_Pipeline();
                OpenGL_Line(deformable_body_collection.particles.X(p1),deformable_body_collection.particles.X(p2));}}
        OpenGL_End();
        glPopAttrib();}

    for(int i=0;i<boundary_surface_objects.m;i++){
        if(!active_list(i)) continue;
        glPushName(i);
        //if(embedded_surface_objects(i) && display_mode==3){glPushName(3);embedded_surface_objects(i)->Display();glPopName();}
        if(boundary_surface_objects(i) && display_boundary_surface_objects){
          boundary_surface_objects(i)->wireframe_only=(display_mode==3);glPushName(5);boundary_surface_objects(i)->Display();glPopName();}
        if(hard_bound_boundary_surface_objects(i) && display_hard_bound_boundary_surface_objects){
            hard_bound_boundary_surface_objects(i)->wireframe_only=(display_mode==3);glPushName(5);hard_bound_boundary_surface_objects(i)->Display();glPopName();}
        glPopName();}
    if(slice && slice->Is_Slice_Mode()) glPopAttrib();

    if(interaction_pair_display_mode){
        glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
        glEnable(GL_BLEND);
        glDisable(GL_LIGHTING);

        // visualize point face interactions
        if(interaction_pair_display_mode==0 || interaction_pair_display_mode==1){
            OpenGL_Begin(GL_LINES);
            for(int k=0;k<point_triangle_interaction_pairs.Size();k++){
                const REPULSION_PAIR<TV>& pair=point_triangle_interaction_pairs(k);
                INDIRECT_ARRAY<const ARRAY_VIEW<TV>,VECTOR<int,4>&> X(deformable_body_collection.particles.X,pair.nodes);
                glColor3f(.5f,1,.5f);
                OpenGL_Line(X(0),PLANE<T>(X(1),X(2),X(3)).Surface(X(0)));
                glColor3f(0,1,1);
                OpenGL_Line(X(0),X(1));OpenGL_Line(X(0),X(2));OpenGL_Line(X(0),X(3));}
            OpenGL_End();
            OpenGL_Begin(GL_TRIANGLES);
            for(int k=0;k<point_triangle_interaction_pairs.Size();k++){
                const REPULSION_PAIR<TV>& pair=point_triangle_interaction_pairs(k);
                INDIRECT_ARRAY<const ARRAY_VIEW<TV>,VECTOR<int,4>&> X(deformable_body_collection.particles.X,pair.nodes);
                glColor4f(0,.6f,.8f,.5f);
                OpenGL_Triangle(X(0),X(2),X(1));
                OpenGL_Triangle(X(0),X(1),X(3));
                OpenGL_Triangle(X(0),X(3),X(2));}
            OpenGL_End();}

        // visualize edge edge interactions
        if(interaction_pair_display_mode==0 || interaction_pair_display_mode==2){
            OpenGL_Begin(GL_LINES);
            for(int k=0;k<edge_edge_interaction_pairs.Size();k++){
                const REPULSION_PAIR<TV>& pair=edge_edge_interaction_pairs(k);
                INDIRECT_ARRAY<const ARRAY_VIEW<TV>,VECTOR<int,4>&> X(deformable_body_collection.particles.X,pair.nodes);
                glColor3f(1,1,0);
                VECTOR<T,2> weights;SEGMENT_3D<T>(X(0),X(1)).Shortest_Vector_Between_Segments(SEGMENT_3D<T>(X(2),X(3)),weights);
                OpenGL_Line((1-weights.x)*X(0)+weights.x*X(1),(1-weights.y)*X(2)+weights.y*X(3));
                //OpenGL_Line(X(0),X(2));OpenGL_Line(X(0),X(3));OpenGL_Line(X(1),X(2));OpenGL_Line(X(1),X(3));
                glColor3f(1,.5,0);
                OpenGL_Line(X(0),X(1));OpenGL_Line(X(2),X(3));}
            OpenGL_End();}
        glPopAttrib();}

    //Draw force data
    if(display_forces_mode){
        glPushAttrib(GL_ENABLE_BIT|GL_CURRENT_BIT);
        glEnable(GL_BLEND);
        glDisable(GL_LIGHTING);
        OpenGL_Begin(GL_LINES);
        for(int k=0;k<force_data_list.Size();k++){
            const FORCE_DATA<TV>& force_data=force_data_list(k);
            if(display_forces_mode==0 && force_data.name!="LINEAR_SPRINGS") continue;
            else if(display_forces_mode==1 && force_data.name!="TRIANGLE_BENDING_SPRINGS") continue;
            else if(display_forces_mode==2 && force_data.name!="LINEAR_ALTITUDE_SPRINGS_3D") continue;
            else if(display_forces_mode==3 && force_data.name!="SCALED_SOLIDS_FORCES") continue;

            OPENGL_COLOR force_color=color_map_forces->Lookup(force_data.state);force_color.Send_To_GL_Pipeline();
            OpenGL_Line(force_data.first_action_point,force_data.second_action_point);}
        OpenGL_End();
        glPopAttrib();}
}
//#####################################################################
// Function Cycle_Display_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Cycle_Display_Mode()
{
    display_mode=(display_mode+1)%4;
    Update_Velocity_Field();
}
//#####################################################################
// Function Cycle_Cutaway_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Cycle_Cutaway_Mode()
{
    for(int i=0;i<tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i))tetrahedralized_volume_objects(i)->Cycle_Cutaway_Mode();
}
//#####################################################################
// Function Show_Only_First
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Show_Only_First()
{
    active_list.Fill(false);
    active_list(0)=true;
    Update_Velocity_Field();
}
//#####################################################################
// Function Decrease_Cutaway_Fraction
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Decrease_Cutaway_Fraction()
{
    for(int i=0;i<tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i))tetrahedralized_volume_objects(i)->Decrease_Cutaway_Fraction();
}
//#####################################################################
// Function Increase_Cutaway_Fraction
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Increase_Cutaway_Fraction()
{
    for(int i=0;i<tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i))tetrahedralized_volume_objects(i)->Increase_Cutaway_Fraction();
}
//#####################################################################
// Function Create_One_Big_Triangulated_Surface_And_Write_To_File
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Create_One_Big_Triangulated_Surface_And_Write_To_File()
{
    GEOMETRY_PARTICLES<VECTOR<T,3> >& particles=deformable_body_collection.particles;TRIANGLE_MESH mesh;
    TRIANGULATED_SURFACE<T> triangulated_surface(mesh,particles);
    for(int index=0;index<deformable_body_collection.structures.m;index++){
        LOG::cout<<"Adding "<<index<<"th triangulated surface out of "<<deformable_body_collection.structures.m<<" surfaces"<<std::endl;
        TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(index));
        if(!tetrahedralized_volume) continue;
        tetrahedralized_volume->Initialize_Triangulated_Surface();
        TRIANGULATED_SURFACE<T>& individual_boundary_triangulated_surface=*tetrahedralized_volume->triangulated_surface;
        int individual_number_of_triangles=individual_boundary_triangulated_surface.mesh.elements.m;
        for(int t=0;t<individual_number_of_triangles;t++) mesh.elements.Append(individual_boundary_triangulated_surface.mesh.elements(t));}
    mesh.number_nodes=particles.Size();
    FILE_UTILITIES::Write_To_File<T>("one_big_tri_surface.tri",triangulated_surface);
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Use_Bounding_Box() const
{
    return draw && valid && deformable_body_collection.structures.m>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Bounding_Box() const
{
    RANGE<VECTOR<T,3> > box=RANGE<VECTOR<T,3> >::Empty_Box();
    if(draw && valid && deformable_body_collection.structures.m>0){
        for(int i=0;i<segmented_curve_objects.m;i++) if(segmented_curve_objects(i))box.Enlarge_To_Include_Box(segmented_curve_objects(i)->Bounding_Box());
        for(int i=0;i<triangulated_surface_objects.m;i++) if(triangulated_surface_objects(i))box.Enlarge_To_Include_Box(triangulated_surface_objects(i)->Bounding_Box());
        for(int i=0;i<tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i))box.Enlarge_To_Include_Box(tetrahedralized_volume_objects(i)->Bounding_Box());
        for(int i=0;i<hexahedralized_volume_objects.m;i++) if(hexahedralized_volume_objects(i))box.Enlarge_To_Include_Box(hexahedralized_volume_objects(i)->Bounding_Box());        
        for(int i=0;i<free_particles_objects.m;i++) if(free_particles_objects(i)) box.Enlarge_To_Include_Box(free_particles_objects(i)->Bounding_Box());
        for(int i=0;i<boundary_surface_objects.m;i++) if(boundary_surface_objects(i))box.Enlarge_To_Include_Box(boundary_surface_objects(i)->Bounding_Box());}
    return box;
}
//#####################################################################
// Function Highlight_Particle
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Highlight_Particle()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Enter Particle Number: ",{[this](){Highlight_Particle_Response();},""});
}
//#####################################################################
// Function Toggle_Active_Value
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Toggle_Active_Value()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Enter Component Number: ",{[this](){Toggle_Active_Value_Response();},""});
}
//#####################################################################
// Function Toggle_Hide_Unselected
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Toggle_Hide_Unselected()
{
    if(real_selection) hide_unselected=!hide_unselected;
}
//#####################################################################
// Function Toggle_Draw_Interior
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Toggle_Draw_Interior()
{
    for(int i=0;i<tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i)) tetrahedralized_volume_objects(i)->Toggle_Boundary_Only();
    for(int i=0;i<hexahedralized_volume_objects.m;i++) if(hexahedralized_volume_objects(i) && active_list(i)) hexahedralized_volume_objects(i)->Toggle_Boundary_Only();
}
//#####################################################################
// Function Toggle_Draw_Subsets
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Toggle_Draw_Subsets()
{
    for(int i=0;i<tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i))
        tetrahedralized_volume_objects(i)->draw_subsets=!tetrahedralized_volume_objects(i)->draw_subsets;
    for(int i=0;i<hexahedralized_volume_objects.m;i++) if(hexahedralized_volume_objects(i) && active_list(i))
        hexahedralized_volume_objects(i)->draw_subsets=!hexahedralized_volume_objects(i)->draw_subsets;
}
//#####################################################################
// Function Toggle_Use_Active_List
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Toggle_Use_Active_List()
{
    if(active_list.m<1) return;
    use_active_list=!use_active_list;if(!use_active_list) active_list.Fill(true);
    else{active_list.Fill(false);incremented_active_object=1;active_list(incremented_active_object)=true;}
    Update_Velocity_Field();
}
//#####################################################################
// Function Select_Segment
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Toggle_Selection_Mode()
{
    if(!real_selection) return;
    assert(real_selection->body_selection);
    if(real_selection->body_selection->type == OPENGL_SELECTION<T>::SEGMENTED_CURVE_3D){
        delete real_selection->body_selection;
        real_selection->body_selection=real_selection->saved_selection;}
    else if(real_selection->body_selection->type == OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_3D){
        real_selection->saved_selection=real_selection->body_selection;
        real_selection->body_selection=((OPENGL_SEGMENTED_CURVE_3D<T>*)real_selection->subobject)->Get_Curve_Selection(
                ((OPENGL_SELECTION_SEGMENTED_CURVE_3D<T>*)real_selection->body_selection)->index);}
    else return;
    Highlight_Selection(real_selection);
    glutPostRedisplay();
}
//#####################################################################
// Function Increment_Active_Object
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Increment_Active_Object()
{
    if(active_list.m<1 || !use_active_list) return;
    active_list(incremented_active_object)=false;incremented_active_object++;
    if(incremented_active_object>active_list.m) incremented_active_object=0;
    active_list(incremented_active_object)=true;
    Update_Velocity_Field();
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>* selection=0;
    if(buffer_size>=1){
        selection=new OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>(this);
        selection->body_index=buffer[0];selection->subobject_type=buffer[1];
        switch(selection->subobject_type){
            case 1:selection->subobject=segmented_curve_objects(buffer[0]);break;
            case 2:selection->subobject=triangulated_surface_objects(buffer[0]);break;
            case 3:selection->subobject=tetrahedralized_volume_objects(buffer[0]);break;
            case 4:selection->subobject=embedded_surface_objects(buffer[0]);break;
            case 5:selection->subobject=boundary_surface_objects(buffer[0]);break;
            case 6:selection->subobject=free_particles_objects(buffer[0]);break;
            default:PHYSBAM_FATAL_ERROR();}
        selection->body_selection=selection->subobject->Get_Selection(&buffer[2],buffer_size-2);
        if(!selection->body_selection){delete selection;selection=0;}}
    return selection;
}
//#####################################################################
// Function Set_Selection
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Set_Selection(OPENGL_SELECTION<T>* selection)
{
    if(selection->type!=OPENGL_SELECTION<T>::COMPONENT_DEFORMABLE_COLLECTION_3D) return;
    real_selection=(OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>*)selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    if(selection->type!=OPENGL_SELECTION<T>::COMPONENT_DEFORMABLE_COLLECTION_3D) return;
    OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>* real_selection=(OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>*)selection;
    if(selection->hide){
        active_list(real_selection->body_index)=false;
        Update_Velocity_Field();}
    else real_selection->subobject->Highlight_Selection(real_selection->body_selection);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Clear_Highlight()
{
    for(int i=0;i<segmented_curve_objects.m;i++) if(segmented_curve_objects(i) && active_list(i)) segmented_curve_objects(i)->Clear_Highlight();
    for(int i=0;i<triangulated_surface_objects.m;i++) if(triangulated_surface_objects(i) && active_list(i)) triangulated_surface_objects(i)->Clear_Highlight();
    for(int i=0;i<tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i)) tetrahedralized_volume_objects(i)->Clear_Highlight();
    for(int i=0;i<hexahedralized_volume_objects.m;i++) if(hexahedralized_volume_objects(i) && active_list(i)) hexahedralized_volume_objects(i)->Clear_Highlight();
    for(int i=0;i<free_particles_objects.m;i++) if(free_particles_objects(i) && active_list(i))free_particles_objects(i)->Clear_Highlight();
    for(int i=0;i<boundary_surface_objects.m;i++) if(boundary_surface_objects(i) && active_list(i))boundary_surface_objects(i)->Clear_Highlight();
    for(int i=0;i<embedded_surface_objects.m;i++) if(embedded_surface_objects(i) && active_list(i)) embedded_surface_objects(i)->Clear_Highlight();
    for(int i=0;i<hard_bound_boundary_surface_objects.m;i++) if(hard_bound_boundary_surface_objects(i) && active_list(i))hard_bound_boundary_surface_objects(i)->Clear_Highlight();
    if(hide_unselected) hide_unselected=false;
    real_selection=0;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION<T>* selection) const
{
    if(selection && selection->type==OPENGL_SELECTION<T>::COMPONENT_DEFORMABLE_COLLECTION_3D && selection->object==this){
        OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>* real_selection=(OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>*)selection;
        output_stream<<"Deformable object "<<real_selection->body_index<<std::endl;
        real_selection->subobject->Print_Selection_Info(output_stream,real_selection->body_selection);}
}
//#####################################################################
// Function Toggle_Active_Value_Response
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Toggle_Active_Value_Response()
{
    bool start_val=true;
    if(!OPENGL_WORLD<T>::Singleton()->prompt_response.empty()){
        int index;
        std::istringstream sstream(OPENGL_WORLD<T>::Singleton()->prompt_response);
        sstream>>index;
        if(index>=0) {
            if(active_list.m<=index) active_list.Resize(index+1,true,true,start_val);
            active_list(index)=!active_list(index);
            Update_Velocity_Field();}}
}
//#####################################################################
// Function Highlight_Particle_Response
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Highlight_Particle_Response()
{
    if(!OPENGL_WORLD<T>::Singleton()->prompt_response.empty()){
        int index=-1;std::istringstream sstream(OPENGL_WORLD<T>::Singleton()->prompt_response);sstream>>index;
        if(index>=0 && index<deformable_body_collection.particles.Size()) selected_vertex=index;}
    Reinitialize(true);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T>
RANGE<VECTOR<T,3> > OPENGL_SELECTION_COMPONENT_DEFORMABLE_COLLECTION_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object && body_selection);
    return object->World_Space_Box(body_selection->Bounding_Box());
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Set_Vector_Size(double size)
{
    velocity_field.size=size;
}
//#####################################################################
// Function Toggle_Velocity_Vectors
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Toggle_Velocity_Vectors()
{
    draw_velocity_vectors=!draw_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Increase_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Decrease_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Update_Velocity_Field()
{
    if(!draw_velocity_vectors) return;
    bool display_triangulated_surface_objects,display_tetrahedralized_volume_objects,display_hexahedralized_volume_objects,
         display_boundary_surface_objects,display_hard_bound_boundary_surface_objects,display_free_particles_objects;
    Set_Display_Modes(display_triangulated_surface_objects,display_tetrahedralized_volume_objects,
            display_hexahedralized_volume_objects,display_boundary_surface_objects,display_hard_bound_boundary_surface_objects,display_free_particles_objects);

    velocity_vectors.Resize(deformable_body_collection.particles.V.m);
    velocity_vectors.Fill(TV());
    for(int i=0;i<segmented_curve_objects.m;i++){
        if(!active_list(i)) continue;
        if(segmented_curve_objects(i))
            velocity_vectors.Subset(segmented_curve_objects(i)->curve.mesh.elements.Flattened())=
                deformable_body_collection.particles.V.Subset(segmented_curve_objects(i)->curve.mesh.elements.Flattened());
        if(triangulated_surface_objects(i) && display_triangulated_surface_objects)
            velocity_vectors.Subset(triangulated_surface_objects(i)->surface.mesh.elements.Flattened())=
                deformable_body_collection.particles.V.Subset(triangulated_surface_objects(i)->surface.mesh.elements.Flattened());
        if(tetrahedralized_volume_objects(i) && display_tetrahedralized_volume_objects)
            velocity_vectors.Subset(tetrahedralized_volume_objects(i)->mesh->elements.Flattened())=
                deformable_body_collection.particles.V.Subset(tetrahedralized_volume_objects(i)->mesh->elements.Flattened());
        if(hexahedralized_volume_objects(i) && display_hexahedralized_volume_objects)
            velocity_vectors.Subset(hexahedralized_volume_objects(i)->hexahedron_mesh->elements.Flattened())=
                deformable_body_collection.particles.V.Subset(hexahedralized_volume_objects(i)->hexahedron_mesh->elements.Flattened());
        if(free_particles_objects(i) && display_free_particles_objects)
            velocity_vectors.Subset(free_particles_objects(i)->points.indices)=
                deformable_body_collection.particles.V.Subset(free_particles_objects(i)->points.indices);}

    positions=deformable_body_collection.particles.X;
}
//#####################################################################
// Function Cycle_Relative_Velocity_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Cycle_Relative_Velocity_Mode()
{
    display_relative_velocity_mode=(display_relative_velocity_mode+1)%number_of_segmented_curve;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION<T>* old_selection,bool& delete_selection)
{
    if(old_selection && old_selection->object==this && invalidate_deformable_objects_selection_each_frame) delete_selection=true;
    return 0;
}
//#####################################################################
// Function Toggle_Differentiate_Inverted
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Toggle_Differentiate_Inverted()
{
    for(int i=0;i<tetrahedralized_volume_objects.m;i++) if(tetrahedralized_volume_objects(i) && active_list(i)) tetrahedralized_volume_objects(i)->Toggle_Differentiate_Inverted();
    for(int i=0;i<hexahedralized_volume_objects.m;i++) if(hexahedralized_volume_objects(i) && active_list(i)) hexahedralized_volume_objects(i)->Toggle_Differentiate_Inverted();
}
//#####################################################################
// Function Create_Hard_Bound_Boundary_Surface
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>& OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Create_Hard_Bound_Boundary_Surface(TRIANGULATED_SURFACE<T>& boundary_surface)
{
    TRIANGULATED_SURFACE<T>& hard_bound_boundary_surface=*TRIANGULATED_SURFACE<T>::Create(boundary_surface.particles);
    ARRAY<int> particle_map(IDENTITY_ARRAY<>(boundary_surface.particles.Size()));
#if 0 // TODO: Fix me
    for(int b=0;b<deformable_body_collection.soft_bindings.bindings.m;b++){VECTOR<int,2>& binding=deformable_body_collection.soft_bindings.bindings(b);particle_map(binding.x)=binding.y;}
#endif
    for(int t=0;t<boundary_surface.mesh.elements.m;t++) hard_bound_boundary_surface.mesh.elements.Append(VECTOR<int,3>::Map(particle_map,boundary_surface.mesh.elements(t)));
    return hard_bound_boundary_surface;
}
//#####################################################################
// Function Cycle_Hard_Bound_Surface_Display_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Cycle_Hard_Bound_Surface_Display_Mode()
{
    display_hard_bound_surface_mode=(display_hard_bound_surface_mode+1)%3;
    display_soft_bound_surface_mode=!has_embedded_objects && has_soft_bindings?display_hard_bound_surface_mode:1; // TODO: fix names
}
//#####################################################################
// Function Cycle_Forces_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Cycle_Forces_Mode()
{
    display_forces_mode=(display_forces_mode+1)%6;
    if(display_forces_mode==0) LOG::cout<<"Displaying no forces"<<std::endl;
    else if(display_forces_mode==1) LOG::cout<<"Displaying forces for LINEAR_SPRINGS"<<std::endl;
    else if(display_forces_mode==2) LOG::cout<<"Displaying forces for TRIANGLE_BENDING_SPRINGS"<<std::endl;
    else if(display_forces_mode==3) LOG::cout<<"Displaying forces for LINEAR_ALTITUDE_SPRINGS_3D"<<std::endl;
    else if(display_forces_mode==4) LOG::cout<<"Displaying forces for SCALED_SOLIDS_FORCES"<<std::endl;
    else LOG::cout<<"Displaying all forces"<<std::endl;
}
//#####################################################################
// Function Cycle_Interaction_Pair_Display_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>::
Cycle_Interaction_Pair_Display_Mode()
{
    if(!interaction_pair_display_mode && !point_triangle_interaction_pairs.m && !edge_edge_interaction_pairs.m){
        std::string file=LOG::sprintf("%s/%d/interaction_pairs",prefix.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(file))
            FILE_UTILITIES::Read_From_File(stream_type,file,point_triangle_interaction_pairs,edge_edge_interaction_pairs);}
    interaction_pair_display_mode=(interaction_pair_display_mode+1)%4;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<float>;
template class OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<double>;
}
