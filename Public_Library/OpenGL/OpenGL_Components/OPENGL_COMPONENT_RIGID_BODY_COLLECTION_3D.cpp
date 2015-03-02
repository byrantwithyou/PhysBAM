//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Joints/ANGLE_JOINT.h>
#include <Rigids/Joints/JOINT_MESH.h>
#include <Rigids/Joints/NORMAL_JOINT.h>
#include <Rigids/Joints/POINT_JOINT.h>
#include <Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
#include <Rigids/Joints/RIGID_JOINT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <OpenGL/OpenGL/OPENGL_AXES.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD.h>
#include <sstream>
#include <string>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(STREAM_TYPE stream_type,const std::string& basedir_input,bool use_display_lists)
    :OPENGL_COMPONENT<T>(stream_type,"Rigid Geometry Collection"),basedir(basedir_input),use_display_lists(use_display_lists),frame_loaded(-1),valid(false),
    rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(0)),articulated_rigid_body(0),
    velocity_field(stream_type,velocity_vectors,positions,OPENGL_COLOR::Cyan(),.25,true,true),
    angular_velocity_field(stream_type,angular_velocity_vectors,positions,OPENGL_COLOR::Magenta(),.25,true,true),need_destroy_rigid_body_collection(true),one_sided(false),
    front_color_map(0),back_color_map(0),
    current_selection(0)
{
    Initialize();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(STREAM_TYPE stream_type,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir_input,bool use_display_lists)
    :OPENGL_COMPONENT<T>(stream_type,"Rigid Geometry Collection"),basedir(basedir_input),use_display_lists(use_display_lists),frame_loaded(-1),valid(false),
    rigid_body_collection(rigid_body_collection),articulated_rigid_body(0),
    velocity_field(stream_type,velocity_vectors,positions,OPENGL_COLOR::Cyan(),.25,true,true),
    angular_velocity_field(stream_type,angular_velocity_vectors,positions,OPENGL_COLOR::Magenta(),.25,true,true),need_destroy_rigid_body_collection(false),one_sided(false),
    front_color_map(0),back_color_map(0),current_selection(0)
{
    Initialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D()
{
    delete front_color_map;
    delete back_color_map;
    opengl_triangulated_surface.Delete_Pointers_And_Clean_Memory();
    opengl_tetrahedralized_volume.Delete_Pointers_And_Clean_Memory();
    opengl_levelset.Delete_Pointers_And_Clean_Memory();
    opengl_octree_levelset_surface.Delete_Pointers_And_Clean_Memory();
    opengl_axes.Delete_Pointers_And_Clean_Memory();
    if(need_destroy_rigid_body_collection) delete &rigid_body_collection;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Initialize()
{
    is_animation=true;
    show_object_names=false;
    output_positions=true;
    draw_velocity_vectors=false;
    draw_angular_velocity_vectors=false;
    draw_individual_axes=false;
    draw_triangulated_surface=true;
    draw_tetrahedralized_volume=true;
    draw_implicit_surface=false;
    read_triangulated_surface=true;
    read_implicit_surface=false;
    read_tetrahedralized_volume=false;
    draw_simplicial_object_particles=false;
    draw_articulation_points=false;
    draw_joint_frames=0;
    draw_forces_and_torques=false;
    has_init_destroy_information=true;

    front_color_map=OPENGL_INDEXED_COLOR_MAP::Rigid_Body_Color_Map();
    back_color_map=OPENGL_INDEXED_COLOR_MAP::Rigid_Body_Color_Map();
}
//#####################################################################
// Function Set_Draw_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Set_Draw_Mode(const int mode)
{
    draw_triangulated_surface=(mode&0x01)!=0;
    draw_tetrahedralized_volume=(mode&0x02)!=0;
    draw_implicit_surface=(mode&0x04)!=0;
}
//#####################################################################
// Function Get_Draw_Mode
//#####################################################################
template<class T> int OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Get_Draw_Mode() const
{
    return 0x04*(int)(draw_implicit_surface)+0x02*(int)(draw_tetrahedralized_volume)+0x01*(int)(draw_triangulated_surface);
}
//#####################################################################
// Function Read_Hints
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Read_Hints(const std::string& filename)
{
    ARRAY<OPENGL_RIGID_BODY_HINTS,int> opengl_hints;
    FILE_UTILITIES::Read_From_File(stream_type,filename,opengl_hints);
    for(int i=0;i<opengl_triangulated_surface.Size();i++) if(opengl_triangulated_surface(i) && i<opengl_hints.Size()){
        opengl_triangulated_surface(i)->Set_Front_Material(opengl_hints(i).material);
        use_object_bounding_box(i)=opengl_hints(i).include_bounding_box;}
}
//#####################################################################
// Function Resize_Structures
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Resize_Structures(const int size)
{
    extra_components.Resize(size);
    for(int i=size;i<opengl_triangulated_surface.m;i++){
        delete opengl_triangulated_surface(i);
        delete opengl_tetrahedralized_volume(i);
        delete opengl_levelset(i);
        delete opengl_octree_levelset_surface(i);
        delete opengl_axes(i);}

    opengl_triangulated_surface.Resize(size);
    opengl_tetrahedralized_volume.Resize(size);
    opengl_levelset.Resize(size);
    opengl_octree_levelset_surface.Resize(size);
    opengl_axes.Resize(size);
    draw_object.Resize(size);
    use_object_bounding_box.Resize(size);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Reinitialize(const bool force,const bool read_geometry)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_particles",basedir.c_str(),frame))) return;

        // TODO: currently reads in all structures, should only read in certain kinds based on read_triangulated_surface,read_implicit_surface,read_tetrahedralized_volume
        rigid_body_collection.Read(stream_type,basedir,frame,&needs_init,&needs_destroy);

        std::string arb_state_file=STRING_UTILITIES::string_sprintf("%s/%d/arb_state",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(arb_state_file)){
            if(!articulated_rigid_body) articulated_rigid_body=new ARTICULATED_RIGID_BODY<TV>(rigid_body_collection); // TODO: read in the actual particles
            articulated_rigid_body->Read(stream_type,basedir,frame);
            Initialize();}
        else{delete articulated_rigid_body;articulated_rigid_body=0;}

        if(FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",basedir.c_str(),frame))){
            Read_Articulated_Information(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",basedir.c_str(),frame));}

        // only enlarge array as we read in more geometry to memory
        int max_number_of_bodies=max(extra_components.Size(),rigid_body_collection.rigid_body_particles.Size());
        Resize_Structures(max_number_of_bodies);

        std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_forces_and_torques",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File(stream_type,filename,forces_and_torques);
        else forces_and_torques.Resize(0);
        if(has_init_destroy_information) for(int i=0;i<needs_destroy.m;i++) Destroy_Geometry(needs_destroy(i));

        std::string rigid_body_colors_file=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_colors",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(rigid_body_colors_file)) FILE_UTILITIES::Read_From_File<T>(rigid_body_colors_file,opengl_colors);
        else{opengl_colors.Resize(max_number_of_bodies);opengl_colors.Fill(OPENGL_COLOR::Cyan());}

        // Initialize bodies which have become active
        if(has_init_destroy_information) for(int i=0;i<needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_body_collection.Is_Active(id));
            Create_Geometry(id);}
        else for(int i=0;i<max_number_of_bodies;i++){if(rigid_body_collection.Is_Active(i)) Create_Geometry(i);} // TODO: can we figure out what bodies need_init

        // Only display real bodies (not ghost bodies)
        if(FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/partition",basedir.c_str(),frame))) {
            ARRAY<int> particles_of_this_partition;
            FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/partition",basedir.c_str(),frame),particles_of_this_partition);
            for(int i=0;i<max_number_of_bodies;i++)
                draw_object(i)=false;
            for(int i=0;i<particles_of_this_partition.Size();i++)
                draw_object(particles_of_this_partition(i))=true;}

        // Update active bodies / remove inactive bodies
        for(int id=0;id<rigid_body_collection.rigid_body_particles.Size();id++){
            if(rigid_body_collection.Is_Active(id)){
                Update_Geometry(id);
                RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(id);
                IMPLICIT_OBJECT<TV>* object_space_implicit_object=body.implicit_object?body.implicit_object->object_space_implicit_object:0;
                if(body.name=="ground" || (object_space_implicit_object && typeid(*object_space_implicit_object)==typeid(ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >)))
                    Set_Object_Material(id,OPENGL_COLOR::Ground_Tan((T)1));
                else{
                    if(one_sided) Set_Object_Material(id,front_color_map->Lookup(Value(id)));
                    else Set_Object_Material(id,front_color_map->Lookup(Value(id)),back_color_map->Lookup(Value(id)));}}
            else Destroy_Geometry(id);}
        for(int id=rigid_body_collection.rigid_body_particles.Size();id<opengl_triangulated_surface.Size();id++) Destroy_Geometry(id);

        frame_loaded=frame;
        valid=true;}

    Update_Object_Labels();

    if(use_display_lists) Initialize_Display_Lists();
}
//#####################################################################
// Function Reinitialize_Without_Files
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Reinitialize_Without_Files(const bool force)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        LOG::cout<<"Reinit being called\n";
        valid=false;

        ARRAY<int> needs_init;
        for(int i=0;i<articulated_rigid_body->rigid_body_collection.rigid_body_particles.Size();i++) if(articulated_rigid_body->rigid_body_collection.Is_Active(i)) needs_init.Append(i);

        // only enlarge array as we read in more geometry to memory
        int max_number_of_bodies=max(extra_components.Size(),rigid_body_collection.rigid_body_particles.Size());
        Resize_Structures(max_number_of_bodies);

        // Initialize bodies which have become active
        for(int i=0;i<needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_body_collection.Is_Active(id));
            Create_Geometry(id);}
        
        // Update active bodies / remove inactive bodies
        for(int id=0;id<rigid_body_collection.rigid_body_particles.Size();id++){
            if(rigid_body_collection.Is_Active(id)) Update_Geometry(id);
            else Destroy_Geometry(id);}
        for(int id=rigid_body_collection.rigid_body_particles.Size();id<opengl_triangulated_surface.Size();id++) Destroy_Geometry(id);

        valid=true;}

    Update_Object_Labels();

    if(use_display_lists) Initialize_Display_Lists();
}
//#####################################################################
// Function Initialize_One_Body
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Initialize_One_Body(const int body_id,const bool force)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        LOG::cout<<"Init for one body being called\n";
        valid=false;

        // only enlarge array as we read in more geometry to memory
        int max_number_of_bodies=max(opengl_triangulated_surface.Size(),extra_components.Size(),rigid_body_collection.rigid_body_particles.Size());
        Resize_Structures(max_number_of_bodies);
        // only enlarge array as we read in more geometry to memory
        opengl_colors.Resize(max_number_of_bodies);opengl_colors.Fill(OPENGL_COLOR::Cyan());
        opengl_triangulated_surface.Resize(max_number_of_bodies);
        opengl_tetrahedralized_volume.Resize(max_number_of_bodies);
        opengl_levelset.Resize(max_number_of_bodies);
        opengl_octree_levelset_surface.Resize(max_number_of_bodies);
        //extra_components.Resize(max_number_of_bodies);
        opengl_axes.Resize(max_number_of_bodies);
        draw_object.Resize(max_number_of_bodies);
        use_object_bounding_box.Resize(max_number_of_bodies);

        // Initialize bodies which have become active
        PHYSBAM_ASSERT(rigid_body_collection.Is_Active(body_id));
        Create_Geometry(body_id);

        // Update active bodies / remove inactive bodies
        for(int id=0;id<rigid_body_collection.rigid_body_particles.Size();id++){
            if(rigid_body_collection.Is_Active(id)) Update_Geometry(id);
            else Destroy_Geometry(id);}
        for(int id=rigid_body_collection.rigid_body_particles.Size();id<opengl_triangulated_surface.Size();id++) Destroy_Geometry(id);

        valid=true;}

    Update_Object_Labels();

    if(use_display_lists) Initialize_Display_Lists();
}
//#####################################################################
// Function Update_Bodies
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Update_Bodies(const bool update_arb_points)
{
    if(articulated_rigid_body && update_arb_points) Update_Articulation_Points();
    // Update active bodies / remove inactive bodies
    for(int id=0;id<rigid_body_collection.rigid_body_particles.Size();id++){
        if(rigid_body_collection.Is_Active(id)) Update_Geometry(id);
        else Destroy_Geometry(id);}
    for(int id=rigid_body_collection.rigid_body_particles.Size();id<opengl_triangulated_surface.Size();id++) Destroy_Geometry(id);
    Update_Object_Labels();
}
//#####################################################################
// Function Create_Geometry
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Create_Geometry(const int id)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(id);
    draw_object(id)=true;use_object_bounding_box(id)=true;
    if(rigid_body.name=="ground") use_object_bounding_box(id)=false; // don't use the ground bounding box
    if(!opengl_axes(id)){opengl_axes(id)=new OPENGL_AXES<T>(stream_type);}

    // add tetrahedralized volume
    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=rigid_body.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>();
    if(tetrahedralized_volume && !opengl_tetrahedralized_volume(id)){
        tetrahedralized_volume->mesh.Initialize_Triangle_Mesh();
        opengl_tetrahedralized_volume(id)=new OPENGL_TETRAHEDRALIZED_VOLUME<T>(stream_type,&tetrahedralized_volume->mesh,&tetrahedralized_volume->particles,
            OPENGL_MATERIAL::Metal(OPENGL_COLOR::Magenta(1,1)),OPENGL_MATERIAL::Metal(OPENGL_COLOR::Yellow(1,1)),true);
        opengl_tetrahedralized_volume(id)->Enslave_Transform_To(*opengl_axes(id));}

    // add implicit object
    if(rigid_body.implicit_object){
        IMPLICIT_OBJECT<TV>* object_space_implicit_object=rigid_body.implicit_object->object_space_implicit_object;
        if(LEVELSET_IMPLICIT_OBJECT<TV>* levelset_implicit_object=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(object_space_implicit_object)){
            if(!opengl_levelset(id)){
                opengl_levelset(id)=new OPENGL_LEVELSET_MULTIVIEW<T>(stream_type);
                opengl_levelset(id)->Set_Levelset(levelset_implicit_object->levelset);
                opengl_levelset(id)->Set_Slice(slice);
                opengl_levelset(id)->Generate_Triangulated_Surface();
                opengl_levelset(id)->Enslave_Transform_To(*opengl_axes(id));}}
        if(IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* implicit_object_transformed=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(object_space_implicit_object)){
            if(!opengl_levelset(id)){
                opengl_levelset(id)=new OPENGL_LEVELSET_MULTIVIEW<T>(stream_type);
                opengl_levelset(id)->Set_Levelset(dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(implicit_object_transformed->object_space_implicit_object)->levelset);
                opengl_levelset(id)->Set_Slice(slice);
                opengl_levelset(id)->Generate_Triangulated_Surface();
                opengl_levelset(id)->Enslave_Transform_To(*opengl_axes(id));
                opengl_levelset(id)->implicit_object_transform=implicit_object_transformed->transform;}}
        else if(!rigid_body.simplicial_object && !opengl_triangulated_surface(id)){
            if(typeid(*object_space_implicit_object)==typeid(ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >))
                use_object_bounding_box(id)=false; // don't use the ground bounding box
            TRIANGULATED_SURFACE<T>* surface=TESSELLATION::Generate_Triangles(*object_space_implicit_object);
            if(surface) rigid_body.Add_Structure(*surface);}}

    // add triangulated surface
    if(rigid_body.simplicial_object && !opengl_triangulated_surface(id)){
        LOG::cout<<"name = "<<rigid_body.name<<std::endl;
        OPENGL_COLOR color=opengl_colors(id);
        opengl_triangulated_surface(id)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,*rigid_body.simplicial_object,false,OPENGL_MATERIAL::Plastic(color),
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Green()));
        opengl_triangulated_surface(id)->Enslave_Transform_To(*opengl_axes(id));
        opengl_triangulated_surface(id)->draw_particles=draw_simplicial_object_particles;}

    // add extra components
    if(opengl_tetrahedralized_volume(id)){
        std::string filename_pattern=STRING_UTILITIES::string_sprintf("%s/accumulated_impulses_%d.%%d",basedir.c_str(),id);
        if(FILE_UTILITIES::Frame_File_Exists(filename_pattern,frame)){
            LOG::cout<<"Adding accumulated impulses to rigid body "<<id<<std::endl;
            OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>* component=
                new OPENGL_COMPONENT_TETRAHEDRALIZED_VOLUME_BASED_VECTOR_FIELD<T>(stream_type,*tetrahedralized_volume,filename_pattern);
            component->opengl_vector_field.Enslave_Transform_To(*opengl_axes(id));
            extra_components(id).Append(component);}}
}
//#####################################################################
// Function Update_Geometry
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Update_Geometry(const int id)
{
    if(opengl_axes(id)) *opengl_axes(id)->frame=FRAME<VECTOR<T,3> >(rigid_body_collection.Rigid_Body(id).Frame());
    if(opengl_tetrahedralized_volume(id)){
        std::string color_map_filename=STRING_UTILITIES::string_sprintf("%s/%d/stress_map_of_tetrahedralized_volume_%d",basedir.c_str(),frame,id);
        if(FILE_UTILITIES::File_Exists(color_map_filename)){
            if(!opengl_tetrahedralized_volume(id)->color_map) opengl_tetrahedralized_volume(id)->color_map=new ARRAY<OPENGL_COLOR>;
            FILE_UTILITIES::Read_From_File(stream_type,color_map_filename,*opengl_tetrahedralized_volume(id)->color_map);}
        else if(opengl_tetrahedralized_volume(id)->color_map){delete opengl_tetrahedralized_volume(id)->color_map;opengl_tetrahedralized_volume(id)->color_map=0;}}
    RIGID_BODY<TV> &rigid_body=rigid_body_collection.Rigid_Body(id);
    if(rigid_body.implicit_object && rigid_body.implicit_object->object_space_implicit_object->update_every_frame)
    {
        if(opengl_levelset(id))
            opengl_levelset(id)->Reset_Surface();
        if(opengl_triangulated_surface(id))
        {
            delete opengl_triangulated_surface(id);
            TRIANGULATED_SURFACE<T>* surface=rigid_body.simplicial_object;
            rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(*rigid_body.implicit_object->object_space_implicit_object));
            delete surface;

            OPENGL_COLOR color=opengl_colors(id);
            opengl_triangulated_surface(id)=new OPENGL_TRIANGULATED_SURFACE<T>(stream_type,*rigid_body.simplicial_object,false,OPENGL_MATERIAL::Plastic(color),
                OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Green()));
            opengl_triangulated_surface(id)->Enslave_Transform_To(*opengl_axes(id));
            opengl_triangulated_surface(id)->draw_particles=draw_simplicial_object_particles;
        }
    }
    for(int i=0;i<extra_components(id).m;i++) extra_components(id)(i)->Set_Frame(frame);
}
//#####################################################################
// Function Destroy_Geometry
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Destroy_Geometry(const int id)
{
    if(draw_object.Size()<id)return; // it's possible that we try to delete geometry that we've never added because of substepping in the simulator
    delete opengl_triangulated_surface(id);opengl_triangulated_surface(id)=0;delete opengl_axes(id);opengl_axes(id)=0;
    if(opengl_tetrahedralized_volume(id)){delete opengl_tetrahedralized_volume(id)->color_map;delete opengl_tetrahedralized_volume(id);opengl_tetrahedralized_volume(id)=0;}
    delete opengl_levelset(id);opengl_levelset(id)=0;delete opengl_axes(id);opengl_axes(id)=0;
    extra_components(id).Delete_Pointers_And_Clean_Memory();
    draw_object(id)=false;
}
//#####################################################################
// Function Initialize_Display_Lists
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Initialize_Display_Lists()
{
#if 0
    cout<<"Initializing display lists..."<<endl;
    int i;
    // First create display lists for original triangulated surface
    for(i=0;i<rigid_body_list.triangulated_surface_list.preallocated.m;i++)
        if(!rigid_body_list.triangulated_surface_list.preallocated(i)){
            opengl_triangulated_surface(i)->Create_Display_List();}
    for(i=0;i<rigid_body_list.triangulated_surface_list.preallocated.m;i++){
        int first_instance=rigid_body_list.triangulated_surface_list.preallocated(i);
        if(first_instance>0){
            int first_instance_id=opengl_triangulated_surface(first_instance)->Get_Display_List_Id();
            opengl_triangulated_surface(i)->Use_Display_List(first_instance_id);}}
#endif
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Update_Object_Labels()
{
    int number_of_drawn_bodies=draw_object.Count_Matches(true);
    if(draw_velocity_vectors || draw_angular_velocity_vectors) positions.Resize(number_of_drawn_bodies);
    if(draw_velocity_vectors) velocity_vectors.Resize(number_of_drawn_bodies);
    if(draw_angular_velocity_vectors) angular_velocity_vectors.Resize(number_of_drawn_bodies);

    int idx=0;
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(draw_object(i)){
        if(draw_velocity_vectors || draw_angular_velocity_vectors){idx++;
            positions(idx)=rigid_body_collection.rigid_body_particles.frame(i).t;
            if(draw_velocity_vectors)
                velocity_vectors(idx)=rigid_body_collection.rigid_body_particles.twist(i).linear;
            if(draw_angular_velocity_vectors){
                //rigid_body_collection.Rigid_Body(i).Update_Angular_Velocity();
                angular_velocity_vectors(idx)=rigid_body_collection.rigid_body_particles.twist(i).angular;}}
        if(opengl_triangulated_surface(i)){
            if(output_positions)
                opengl_triangulated_surface(i)->Set_Name(STRING_UTILITIES::string_sprintf("%s <%.3f %.3f %.3f>",rigid_body_collection.Rigid_Body(i).name.c_str(),rigid_body_collection.rigid_body_particles.frame(i).t.x,rigid_body_collection.rigid_body_particles.frame(i).t.y,rigid_body_collection.rigid_body_particles.frame(i).t.z));
            else opengl_triangulated_surface(i)->Set_Name(rigid_body_collection.Rigid_Body(i).name);}}
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_particles",basedir.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    if(draw_input) draw_object.Fill(true);
}
//#####################################################################
// Function Draw_All_Objects
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Draw_All_Objects()
{
    Set_Draw(true);
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Use_Bounding_Box() const
{
    int num_drawable_objects=0;
    for(int i=0;i<opengl_triangulated_surface.Size();i++)
        if(draw_object(i) && use_object_bounding_box(i))
            num_drawable_objects++;
    return draw && num_drawable_objects>0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Bounding_Box() const
{
    RANGE<VECTOR<T,3> > box;
    if(draw)
        for(int i=0;i<opengl_triangulated_surface.Size();i++)
            if(rigid_body_collection.Is_Active(i) && rigid_body_collection.Rigid_Body(i).name!="ground")
                if(draw_object(i) && use_object_bounding_box(i) && opengl_triangulated_surface(i))
                    box=RANGE<VECTOR<T,3> >::Combine(box,opengl_triangulated_surface(i)->Bounding_Box());
    return box;
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION<T>* OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION<T>* selection=0;
    if(buffer_size>=2){
        int body_id(buffer[1]);
        OPENGL_SELECTION<T>* body_selection=0;
        if(buffer[0]==1){ // segmented curve
            PHYSBAM_ASSERT(opengl_triangulated_surface(body_id));
            body_selection=opengl_triangulated_surface(body_id)->Get_Selection(&buffer[2],buffer_size-2);}
        else if(buffer[0]==2){ // tetrahedralized_volume
            PHYSBAM_ASSERT(opengl_tetrahedralized_volume(body_id));
            body_selection=opengl_tetrahedralized_volume(body_id)->Get_Selection(&buffer[2],buffer_size-2);}
        else if(buffer[0]==4){ // articulation joints
            selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>(this,buffer[1]);}
        if(body_selection) selection=new OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D<T>(this,body_id,body_selection);}
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Highlight_Selection(OPENGL_SELECTION<T>* selection)
{
    if(selection->type==OPENGL_SELECTION<T>::COMPONENT_RIGID_BODIES_3D){
        OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D<T>*)selection;
        if(selection->hide) draw_object(real_selection->body_id)=false;
        else if(real_selection->body_selection->type==OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_VERTEX ||
           real_selection->body_selection->type==OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_SEGMENT ||
           real_selection->body_selection->type==OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_TRIANGLE){ // triangulated surface
            if(opengl_triangulated_surface(real_selection->body_id)) // might have become inactive
                opengl_triangulated_surface(real_selection->body_id)->Highlight_Selection(real_selection->body_selection);}
        else if(real_selection->body_selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_VERTEX ||
                real_selection->body_selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_TETRAHEDRON){ // tetrahedralized volume
            if(opengl_tetrahedralized_volume(real_selection->body_id)) // might have become inactive
                opengl_tetrahedralized_volume(real_selection->body_id)->Highlight_Selection(real_selection->body_selection);}}
    delete current_selection;
    current_selection=0;
    if(selection->type==OPENGL_SELECTION<T>::ARTICULATED_RIGID_BODIES_JOINT_3D){
        current_selection=new OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>(this,((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>*)selection)->joint_id);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Clear_Highlight()
{
    for(int i=0;i<opengl_triangulated_surface.Size();i++)
        if(opengl_triangulated_surface(i)) opengl_triangulated_surface(i)->Clear_Highlight();
    for(int i=0;i<opengl_tetrahedralized_volume.Size();i++)
        if(opengl_tetrahedralized_volume(i)) opengl_tetrahedralized_volume(i)->Clear_Highlight();
    delete current_selection;
    current_selection=0;
}
//#####################################################################
// Function Turn_Smooth_Shading_On
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Turn_Smooth_Shading_On()
{
    for(int i=0;i<opengl_triangulated_surface.Size();i++)
        if(opengl_triangulated_surface(i)) opengl_triangulated_surface(i)->Turn_Smooth_Shading_On();
}
//#####################################################################
// Function Turn_Smooth_Shading_Off
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Turn_Smooth_Shading_Off()
{
    for(int i=0;i<opengl_triangulated_surface.Size();i++)
        if(opengl_triangulated_surface(i)) opengl_triangulated_surface(i)->Turn_Smooth_Shading_Off();
}
//#####################################################################
// Function Slice_Has_Changed
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Slice_Has_Changed()
{
    for(int i=0;i<opengl_levelset.Size();i++){
        if(opengl_levelset(i)) opengl_levelset(i)->Set_Slice(slice);
        if(opengl_levelset(i)) opengl_levelset(i)->Slice_Has_Changed();}
}
//#####################################################################
// Function Set_Draw_Object
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Set_Draw_Object(int i,bool draw_it)
{
    draw_object(i)=draw_it;
}
//#####################################################################
// Function Get_Draw_Object
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Get_Draw_Object(int i) const
{
    return draw_object(i);
}
//#####################################################################
// Function Set_Object_Material
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Set_Object_Material(int i,const OPENGL_MATERIAL &front_material_input)
{
    if(!opengl_triangulated_surface(i)) return;
    opengl_triangulated_surface(i)->Set_Front_Material(front_material_input);
    opengl_triangulated_surface(i)->Set_Two_Sided(false);
}
//#####################################################################
// Function Set_Object_Material
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Set_Object_Material(int i,const OPENGL_MATERIAL &front_material_input,const OPENGL_MATERIAL &back_material_input)
{
    if(!opengl_triangulated_surface(i)) return;
    opengl_triangulated_surface(i)->Set_Front_Material(front_material_input);
    opengl_triangulated_surface(i)->Set_Back_Material(back_material_input);
    opengl_triangulated_surface(i)->Set_Two_Sided(true);
}
//#####################################################################
// Function Set_Use_Bounding_Box
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Set_Use_Object_Bounding_Box(int i,bool use_it)
{
    use_object_bounding_box(i)=use_it;
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Set_Vector_Size(double size)
{
    velocity_field.size=size;
    angular_velocity_field.size=size;
}
//#####################################################################
// Function Toggle_Velocity_Vectors
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Velocity_Vectors()
{
    draw_velocity_vectors=!draw_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Angular_Velocity_Vectors
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Angular_Velocity_Vectors()
{
    draw_angular_velocity_vectors=!draw_angular_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Individual_Axes
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Individual_Axes()
{
    draw_individual_axes=!draw_individual_axes;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Output_Positions
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Output_Positions()
{
    output_positions=!output_positions;
    Update_Object_Labels();
}
//#####################################################################
// Function Toggle_Show_Object_Names
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Show_Object_Names()
{
    show_object_names=!show_object_names;
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Increase_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1.1);
    angular_velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Decrease_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1/1.1);
    angular_velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Draw_Values
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Draw_Values()
{
    velocity_field.draw_value=!velocity_field.draw_value;
    angular_velocity_field.draw_value=!angular_velocity_field.draw_value;
}
//#####################################################################
// Function Toggle_Draw_Values
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_One_Sided()
{
    one_sided=!one_sided;
    Reinitialize(true);
}
//#####################################################################
// Function Toggle_Draw_Particles
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Draw_Particles()
{
    draw_simplicial_object_particles=!draw_simplicial_object_particles;
    for(int i=0;i<opengl_triangulated_surface.m;i++)
        if(opengl_triangulated_surface(i)) opengl_triangulated_surface(i)->draw_particles=draw_simplicial_object_particles;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Draw_Mode()
{
    int mode=Get_Draw_Mode();
    mode=(mode+1)%8;
    Set_Draw_Mode(mode);
}
//#####################################################################
// Function Turn_Off_Individual_Smooth_Shading_Prompt
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Turn_Off_Individual_Smooth_Shading_Prompt()
{
    if(!OPENGL_WORLD<T>::Singleton()->prompt_response.empty()){
        int object_id;
        STRING_UTILITIES::String_To_Value(OPENGL_WORLD<T>::Singleton()->prompt_response,object_id);
        if((unsigned)object_id<(unsigned)rigid_body_collection.rigid_body_particles.Size() && opengl_triangulated_surface(object_id))
            opengl_triangulated_surface(object_id)->Turn_Smooth_Shading_Off();}
}
//#####################################################################
// Function Turn_Off_Individual_Smooth_Shading
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Turn_Off_Individual_Smooth_Shading()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Turn off smooth shading for object: ",Turn_Off_Individual_Smooth_Shading_Prompt_CB(),"");
}
//#####################################################################
// Function Manipulate_Individual_Body_Prompt
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Manipulate_Individual_Body_Prompt()
{
    if(!OPENGL_WORLD<T>::Singleton()->prompt_response.empty()){
        int object_id;
        std::string command;
        std::istringstream sstream(OPENGL_WORLD<T>::Singleton()->prompt_response);
        sstream>>command;
        sstream>>object_id;
        if((unsigned)object_id<(unsigned)rigid_body_collection.rigid_body_particles.Size() && opengl_triangulated_surface(object_id)){
            if(command=="s"){
                VECTOR<T,3> scale;
                sstream>>scale;
                LOG::cout<<"Scaling object "<<object_id<<" by "<<scale<<std::endl;
                opengl_triangulated_surface(object_id)->Rescale(scale.x,scale.y,scale.z);}}}
}
//#####################################################################
// Function Manipulate_Individual_Body
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Manipulate_Individual_Body()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Manipulate: ",Manipulate_Individual_Body_Prompt_CB(),"");
}
//#####################################################################
// Selection Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object && body_selection);
    return object->World_Space_Box(body_selection->Bounding_Box());
}
//#####################################################################
// Function Read_Articulated_Information
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Read_Articulated_Information(const std::string& filename)
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    TYPED_ISTREAM typed_input(*input,stream_type);
    // this will need to be changed to reflect multiple articulation points per rigid body
    int numpoints=0;Read_Binary(typed_input,numpoints);articulation_points.Exact_Resize(numpoints);joint_frames.Exact_Resize(numpoints);
    for(int i=0;i<numpoints;i++) Read_Binary(typed_input,articulation_points(i),joint_frames(i));
    Read_Binary(typed_input,projected_COM);delete input;
}
//#####################################################################
// Function Update_Articulation_Points
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Update_Articulation_Points()
{
    int num_points=2*articulated_rigid_body->joint_mesh.Num_Joints();
    articulation_points.Exact_Resize(num_points);
    joint_frames.Exact_Resize(num_points);
    for(int i=0;i<num_points;i+=2){
        int index=i/2;
        JOINT<TV>* joint=articulated_rigid_body->joint_mesh.Joints(index);
        RIGID_BODY<TV>* parent=articulated_rigid_body->Parent(joint->id_number),*child=articulated_rigid_body->Child(joint->id_number);
        articulation_points(i)=parent->World_Space_Point(joint->F_pj().t);articulation_points(i+1)=child->World_Space_Point(joint->F_cj().t);
        joint_frames(i)=parent->Frame()*joint->F_pj();joint_frames(i+1)=child->Frame()*joint->F_cj();}
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Display() const
{
    if(!draw) return;
    GLint mode=0;
    
    glGetIntegerv(GL_RENDER_MODE,&mode);

    if(slice && slice->Is_Slice_Mode()){
        glPushAttrib(GL_ENABLE_BIT);
        slice->Enable_Clip_Planes();}
    if(draw_triangulated_surface){
        glPushName(1);
        for(int i=0;i<opengl_triangulated_surface.Size();i++) if(draw_object(i) && opengl_triangulated_surface(i)){
            glPushName(Value(i));opengl_triangulated_surface(i)->Display();glPopName();}
        glPopName();}
    if(draw_tetrahedralized_volume){
        glPushName(2);
        for(int i=0;i<opengl_tetrahedralized_volume.Size();i++) if(draw_object(i) && opengl_tetrahedralized_volume(i)){
            glPushName(Value(i));opengl_tetrahedralized_volume(i)->Display();glPopName();}
        glPopName();}
    if(draw_implicit_surface){
        glPushName(3);
        int levelset_count=0;
        for(int i=0;i<opengl_levelset.Size();i++) if(draw_object(i)){
            glPushName(Value(i));
            if(opengl_levelset(i)){
                if(++levelset_count>50){
                    PHYSBAM_WARNING("Refusing to draw more than 10 levelsets to save memory.");
                    OPENGL_WORLD<T>::Singleton()->Add_String("WARNING: Refusing to draw more than 10 levelsets to save memory.");
                    break;}
                opengl_levelset(i)->Update();
                opengl_levelset(i)->Display();}
            if(opengl_octree_levelset_surface(i)) opengl_octree_levelset_surface(i)->Display();
            glPopName();}
        glPopName();}
    if(draw_individual_axes)
        for(int i=0;i<opengl_axes.Size();i++){
            if(draw_object(i) && opengl_axes(i)){
                opengl_axes(i)->box.max_corner.x=opengl_axes(i)->box.max_corner.y=opengl_axes(i)->box.max_corner.z=2*rigid_body_collection.Rigid_Body(i).Object_Space_Bounding_Box().Edge_Lengths().Min();
                opengl_axes(i)->Display();}}
    if(draw_velocity_vectors) velocity_field.Display();
    if(draw_angular_velocity_vectors) angular_velocity_field.Display();

    if(show_object_names){
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glColor3f(1,1,1);
        for(int i=0;i<opengl_triangulated_surface.Size();i++)
            if(draw_object(i) && rigid_body_collection.Rigid_Body(i).name.length())
                OpenGL_String(rigid_body_collection.rigid_body_particles.frame(i).t,rigid_body_collection.Rigid_Body(i).name);
        glPopAttrib();}

    // Articulated rigid bodies
    if(articulated_rigid_body){
        if(draw_articulation_points){
            OPENGL_COLOR articulation_point_color(0.8,0.8,0.2),segment_color(0,0,1),com_color(1,0,0);
            if(mode==GL_SELECT){glPushName(4);glPushAttrib(GL_POINT_BIT);glPointSize(OPENGL_PREFERENCES::selection_point_size);}
            for(int i=0;i<articulation_points.m;i++){
                glPushName(i);
                OPENGL_SHAPES::Draw_Dot(articulation_points(i),articulation_point_color,5);
                glPopName();}
            OPENGL_SHAPES::Draw_Dot(projected_COM,com_color,10);
            if(mode!=GL_SELECT) for(int i=0;i<articulation_points.m;i+=2){
                OPENGL_SHAPES::Draw_Segment(articulation_points(i),articulation_points(i+1),segment_color,5);}
            if(mode==GL_SELECT){glPopName();glPopAttrib();}
            if(mode!=GL_SELECT && current_selection && current_selection->type==OPENGL_SELECTION<T>::ARTICULATED_RIGID_BODIES_JOINT_3D){
                int joint_id=((OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>*)current_selection)->joint_id;
                OPENGL_SELECTION<T>::Draw_Highlighted_Vertex(articulation_points(joint_id));}}}

    RANGE<TV> axes_box(RANGE<TV>::Unit_Box()*2);
    //RANGE<TV> axes_box(0,velocity_field.size,0,velocity_field.size,0,velocity_field.size);
    if(draw_joint_frames==1) for(int i=0;i<joint_frames.m;i++)(OPENGL_AXES<T>(stream_type,joint_frames(i),axes_box)).Display();
    else if(draw_joint_frames==2) for(int i=1;i<joint_frames.m;i+=2)(OPENGL_AXES<T>(stream_type,joint_frames(i),axes_box)).Display();
    else if(draw_joint_frames==3) for(int i=0;i<joint_frames.m;i+=2)(OPENGL_AXES<T>(stream_type,joint_frames(i),axes_box)).Display();

    if(draw_forces_and_torques && forces_and_torques.Size()==rigid_body_collection.rigid_body_particles.Size()){
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING);
        T scale=(T)velocity_field.size/24;
        OPENGL_COLOR::Yellow().Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        for(int i=0;i<forces_and_torques.Size();i++) if(rigid_body_collection.Is_Active(i)){
            OpenGL_Line(rigid_body_collection.rigid_body_particles.frame(i).t,rigid_body_collection.rigid_body_particles.frame(i).t+scale*forces_and_torques(i).x);}
        OpenGL_End();
        for(int i=0;i<forces_and_torques.Size();i++) if(rigid_body_collection.Is_Active(i)){
            std::string label=STRING_UTILITIES::string_sprintf("F=%.3f %.3f %.3f, T=%.3f %.3f %.3f",forces_and_torques(i).x.x,forces_and_torques(i).x.y,forces_and_torques(i).x.z,forces_and_torques(i).y.x,forces_and_torques(i).y.y,forces_and_torques(i).y.z);
            OpenGL_String(rigid_body_collection.rigid_body_particles.frame(i).t+scale*forces_and_torques(i).x,label);}
        glPopAttrib();}

    if(slice && slice->Is_Slice_Mode()) glPopAttrib();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION<T>* selection) const
{
    if(!selection || selection->object!=this) return;

    if(selection->type==OPENGL_SELECTION<T>::COMPONENT_RIGID_BODIES_3D){
        OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D<T> *real_selection=(OPENGL_SELECTION_COMPONENT_RIGID_BODY_COLLECTION_3D<T>*)selection;

        output_stream<<"Rigid body "<<real_selection->body_id<<std::endl;

        // handle case of body no longer being active
        if(!rigid_body_collection.Is_Active(real_selection->body_id) || !opengl_triangulated_surface(real_selection->body_id)){
            output_stream<<"INACTIVE"<<std::endl;
            return;}

        RIGID_BODY<TV> *body=const_cast<RIGID_BODY<TV>*>(&rigid_body_collection.Rigid_Body(real_selection->body_id));
        body->Update_Angular_Velocity();

        if(!body->name.empty()){output_stream<<"Name ="<<body->name<<std::endl;}
        output_stream<<"Mass = "<<body->Mass()<<std::endl;
        output_stream<<"Inertia tensor = "<<body->Inertia_Tensor()<<std::endl;
        output_stream<<"Angular momentum = "<<body->Angular_Momentum()<<std::endl;
        output_stream<<"Kinetic energy = "<<body->Kinetic_Energy()<<std::endl;
        output_stream<<std::endl;

        rigid_body_collection.rigid_body_particles.Print(output_stream,real_selection->body_id);

        MATRIX<T,4> body_transform=body->Frame().Matrix();

        if(real_selection->body_selection->type==OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_VERTEX ||
            real_selection->body_selection->type==OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_SEGMENT ||
            real_selection->body_selection->type==OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_TRIANGLE){ // triangulated surface
            if(opengl_triangulated_surface(real_selection->body_id))
                opengl_triangulated_surface(real_selection->body_id)->Print_Selection_Info(output_stream,real_selection->body_selection,&body_transform);
            if(real_selection->body_selection->type==OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_VERTEX){
                VECTOR<T,3> position=body->simplicial_object->particles.X(((OPENGL_SELECTION_TRIANGULATED_SURFACE_VERTEX<T>*)real_selection->body_selection)->index);
                output_stream<<"Pointwise velocity = "<<body->Pointwise_Object_Velocity(body->World_Space_Point(position))<<std::endl;}}
        else if(real_selection->body_selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_VERTEX ||
                real_selection->body_selection->type==OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_TETRAHEDRON){
            if(opengl_tetrahedralized_volume(real_selection->body_id))
                opengl_tetrahedralized_volume(real_selection->body_id)->Print_Selection_Info(output_stream,real_selection->body_selection,&body_transform);}}
    else if(selection->type==OPENGL_SELECTION<T>::ARTICULATED_RIGID_BODIES_JOINT_3D){
        OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T> *real_selection=(OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>*)selection;
        int articulation_id=real_selection->joint_id;int joint_index=(articulation_id+1/2);
        JOINT<TV>& joint=*articulated_rigid_body->joint_mesh.Joints(joint_index);JOINT_ID joint_id=joint.id_number;
        const RIGID_BODY<TV> *parent=articulated_rigid_body->Parent(joint_id),*child=articulated_rigid_body->Child(joint_id);
        output_stream<<"Joint id "<<joint_id<<" ("<<(!joint.name.empty()?joint.name:"UNNAMED")<<")"<<std::endl;
        const char* joint_type_string;
        const std::type_info& type=typeid(joint);
        if(type==typeid(POINT_JOINT<TV>)) joint_type_string="POINT_JOINT";
        else if(type==typeid(RIGID_JOINT<TV>)) joint_type_string="RIGID_JOINT";
        else if(type==typeid(ANGLE_JOINT<TV>)) joint_type_string="ANGLE_JOINT";
        else if(type==typeid(PRISMATIC_TWIST_JOINT<TV>)) joint_type_string="PRISMATIC_TWIST_JOINT";
        else if(type==typeid(NORMAL_JOINT<TV>)) joint_type_string="NORMAL_JOINT";
        else joint_type_string="unknown";
        output_stream<<"Type = "<<joint_type_string<<std::endl;
        output_stream<<"Frame: "<<(articulation_id%2==0?"child":"parent")<<std::endl;
        output_stream<<"Parent = "<<parent->name<<std::endl;
        output_stream<<"Child = "<<child->name<<std::endl;
        output_stream<<"Articulation point "<<articulation_points(articulation_id)<<std::endl;

        FRAME<TV> parent_frame=joint_frames(2*Value(joint_index)-1),child_frame=joint_frames(2*Value(joint_index));
        FRAME<TV> joint_frame=parent_frame.Inverse()*child_frame;
        output_stream<<"Joint translation: "<<joint_frame.t<<std::endl;
        output_stream<<"Joint rotation vector: "<<joint_frame.r.Rotation_Vector()<<std::endl;
        T twist,phi,theta;joint_frame.r.Euler_Angles(twist,phi,theta);
        output_stream<<"Joint Euler angles: twist="<<twist<<", phi="<<phi<<", theta="<<theta<<std::endl;
        TV ap1=parent->World_Space_Point(joint.F_pj().t),ap2=child->World_Space_Point(joint.F_cj().t),location=(T).5*(ap1+ap2);

        TV current_relative_velocity=-RIGID_BODY<TV>::Relative_Velocity(*parent,*child,location); // child w.r.t. parent!
        TV current_relative_angular_velocity=-RIGID_BODY<TV>::Relative_Angular_Velocity(*parent,*child); // child w.r.t. parent!
        output_stream<<"Relative velocity at joint = "<<current_relative_velocity<<std::endl;
        output_stream<<"Relative angular velocity = "<<current_relative_angular_velocity<<std::endl;}
}
//#####################################################################
// Function Toggle_Articulation_Points
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Articulation_Points()
{
    draw_articulation_points=!draw_articulation_points;
}
//#####################################################################
// Function Toggle_Joint_Frames
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Joint_Frames()
{
    draw_joint_frames++;
    draw_joint_frames%=4;
}
//#####################################################################
// Function Toggle_Forces_And_Torques
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Toggle_Forces_And_Torques()
{
    draw_forces_and_torques=!draw_forces_and_torques;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SELECTION_ARTICULATED_RIGID_BODIES_JOINT_3D<T>::
Bounding_Box() const
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<float>;
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<double>;
}

