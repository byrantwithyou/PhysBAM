//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
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
    :OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D(stream_type,*new RIGID_BODY_COLLECTION<TV>(0),basedir_input,use_display_lists)
{
    need_destroy_rigid_body_collection=true;
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
    front_color_map(0),back_color_map(0),selected_joint_id(-1),selected_surface(-1),selected_volume(-1)
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
    opengl_axes.Delete_Pointers_And_Clean_Memory();
    if(need_destroy_rigid_body_collection) delete &rigid_body_collection;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Initialize()
{
    viewer_callbacks.Set("toggle_velocity_vectors",{[this](){Toggle_Velocity_Vectors();},"Toggle velocity vectors"});
    viewer_callbacks.Set("toggle_angular_velocity_vectors",{[this](){Toggle_Angular_Velocity_Vectors();},"Toggle angular velocity vectors"});
    viewer_callbacks.Set("toggle_individual_axes",{[this](){Toggle_Individual_Axes();},"Toggle individual axes"});
    viewer_callbacks.Set("toggle_output_positions",{[this](){Toggle_Output_Positions();},"Toggle output positions"});
    viewer_callbacks.Set("toggle_show_object_names",{[this](){Toggle_Show_Object_Names();},"Toggle show object names"});
    viewer_callbacks.Set("turn_off_individual_smooth_shading",{[this](){Turn_Off_Individual_Smooth_Shading();},"Turn off individual smooth shading"});
    viewer_callbacks.Set("manipulate_individual_body",{[this](){Manipulate_Individual_Body();},"Manipulate individual body"});
    viewer_callbacks.Set("toggle_draw_mode",{[this](){Toggle_Draw_Mode();},"Toggle draw mode"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_draw_values",{[this](){Toggle_Draw_Values();},"Toggle draw values"});
    viewer_callbacks.Set("toggle_one_sided",{[this](){Toggle_One_Sided();},"Toggle one/two sided drawing"});
    viewer_callbacks.Set("toggle_draw_particles",{[this](){Toggle_Draw_Particles();},"Toggle drawing of simplicial object's vertices"});
    viewer_callbacks.Set("toggle_articulation_points",{[this](){Toggle_Articulation_Points();},"Toggle articulation points"});
    viewer_callbacks.Set("toggle_joint_frames",{[this](){Toggle_Joint_Frames();},"Toggle joint frames"});
    viewer_callbacks.Set("toggle_forces_and_torques",{[this](){Toggle_Forces_And_Torques();},"Toggle forces and torques"});

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
        delete opengl_axes(i);}

    opengl_triangulated_surface.Resize(size);
    opengl_tetrahedralized_volume.Resize(size);
    opengl_levelset.Resize(size);
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
        if(!File_Exists(LOG::sprintf("%s/%d/rigid_body_particles",basedir.c_str(),frame))) return;

        // TODO: currently reads in all structures, should only read in certain kinds based on read_triangulated_surface,read_implicit_surface,read_tetrahedralized_volume
        rigid_body_collection.Read(stream_type,basedir,frame,&needs_init,&needs_destroy);

        std::string arb_state_file=LOG::sprintf("%s/%d/arb_state",basedir.c_str(),frame);
        if(File_Exists(arb_state_file)){
            if(!articulated_rigid_body) articulated_rigid_body=new ARTICULATED_RIGID_BODY<TV>(rigid_body_collection); // TODO: read in the actual particles
            articulated_rigid_body->Read(stream_type,basedir,frame);
            Initialize();}
        else{delete articulated_rigid_body;articulated_rigid_body=0;}

        if(File_Exists(LOG::sprintf("%s/%d/arb_info",basedir.c_str(),frame))){
            Read_Articulated_Information(LOG::sprintf("%s/%d/arb_info",basedir.c_str(),frame));}

        // only enlarge array as we read in more geometry to memory
        int max_number_of_bodies=max(extra_components.Size(),rigid_body_collection.rigid_body_particles.Size());
        Resize_Structures(max_number_of_bodies);

        std::string filename=LOG::sprintf("%s/%d/rigid_body_forces_and_torques",basedir.c_str(),frame);
        if(File_Exists(filename)) Read_From_File(stream_type,filename,forces_and_torques);
        else forces_and_torques.Resize(0);
        if(has_init_destroy_information) for(int i=0;i<needs_destroy.m;i++) Destroy_Geometry(needs_destroy(i));

        std::string rigid_body_colors_file=LOG::sprintf("%s/%d/rigid_body_colors",basedir.c_str(),frame);
        if(File_Exists(rigid_body_colors_file)) Read_From_File<T>(rigid_body_colors_file,opengl_colors);
        else{opengl_colors.Resize(max_number_of_bodies);opengl_colors.Fill(OPENGL_COLOR::Cyan());}

        // Initialize bodies which have become active
        if(has_init_destroy_information) for(int i=0;i<needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_body_collection.Is_Active(id));
            Create_Geometry(id);}
        else for(int i=0;i<max_number_of_bodies;i++){if(rigid_body_collection.Is_Active(i)) Create_Geometry(i);} // TODO: can we figure out what bodies need_init

        // Only display real bodies (not ghost bodies)
        if(File_Exists(LOG::sprintf("%s/%d/partition",basedir.c_str(),frame))) {
            ARRAY<int> particles_of_this_partition;
            Read_From_File(stream_type,LOG::sprintf("%s/%d/partition",basedir.c_str(),frame),particles_of_this_partition);
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
        std::string filename_pattern=LOG::sprintf("%s/accumulated_impulses_%d.%%d",basedir.c_str(),id);
        if(Frame_File_Exists(filename_pattern,frame)){
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
        std::string color_map_filename=LOG::sprintf("%s/%d/stress_map_of_tetrahedralized_volume_%d",basedir.c_str(),frame,id);
        if(File_Exists(color_map_filename)){
            if(!opengl_tetrahedralized_volume(id)->color_map) opengl_tetrahedralized_volume(id)->color_map=new ARRAY<OPENGL_COLOR>;
            Read_From_File(stream_type,color_map_filename,*opengl_tetrahedralized_volume(id)->color_map);}
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
                opengl_triangulated_surface(i)->Set_Name(LOG::sprintf("%s <%.3f %.3f %.3f>",rigid_body_collection.Rigid_Body(i).name.c_str(),rigid_body_collection.rigid_body_particles.frame(i).t.x,rigid_body_collection.rigid_body_particles.frame(i).t.y,rigid_body_collection.rigid_body_particles.frame(i).t.z));
            else opengl_triangulated_surface(i)->Set_Name(rigid_body_collection.Rigid_Body(i).name);}}
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Valid_Frame(int frame_input) const
{
    return File_Exists(LOG::sprintf("%s/%d/rigid_body_particles",basedir.c_str(),frame_input));
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
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    const static int priority[]={0,0,10,95};
    if(indices(0)==0) // segmented curve
        return opengl_triangulated_surface(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
    if(indices(0)==1) // tetrahedralized_volume
        return opengl_tetrahedralized_volume(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
    return priority[indices(0)];
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(indices(0)<=1 && (modifiers&GLUT_ACTIVE_CTRL)){
        draw_object(indices(1))=false;
        return false;}
    if(indices(0)==0){
        selected_surface=indices(1);
        return opengl_triangulated_surface(indices(1))->Set_Selection(indices.Array_View(2,indices.m-2),modifiers);}
    else if(indices(0)==1){
        selected_volume=indices(1);
        return opengl_tetrahedralized_volume(indices(1))->Set_Selection(indices.Array_View(2,indices.m-2),modifiers);}
    else if(indices(0)==3) // articulation joints
        selected_joint_id=indices(1);
    else return false;
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Clear_Selection()
{
    if(selected_surface>=0){
        opengl_triangulated_surface(selected_surface)->Clear_Selection();
        selected_surface=-1;}
    if(selected_volume>=0){
        opengl_tetrahedralized_volume(selected_volume)->Clear_Selection();
        selected_volume=-1;}
    selected_joint_id=-1;
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
        String_To_Value(OPENGL_WORLD<T>::Singleton()->prompt_response,object_id);
        if((unsigned)object_id<(unsigned)rigid_body_collection.rigid_body_particles.Size() && opengl_triangulated_surface(object_id))
            opengl_triangulated_surface(object_id)->Turn_Smooth_Shading_Off();}
}
//#####################################################################
// Function Turn_Off_Individual_Smooth_Shading
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Turn_Off_Individual_Smooth_Shading()
{
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Turn off smooth shading for object: ",{[this](){Turn_Off_Individual_Smooth_Shading_Prompt();},""});
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
    OPENGL_WORLD<T>::Singleton()->Prompt_User("Manipulate: ",{[this](){Manipulate_Individual_Body_Prompt();},""});
}
//#####################################################################
// Selection Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Selection_Bounding_Box() const
{
    if(selected_surface>=0) return opengl_triangulated_surface(selected_surface)->Selection_Bounding_Box();
    if(selected_volume>=0) return opengl_tetrahedralized_volume(selected_volume)->Selection_Bounding_Box();
    return RANGE<TV>::Centered_Box();
}
//#####################################################################
// Function Read_Articulated_Information
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Read_Articulated_Information(const std::string& filename)
{
    std::istream* input=Safe_Open_Input(filename);
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
        glPushName(0);
        glPushName(0);
        for(int i=0;i<opengl_triangulated_surface.Size();i++) if(draw_object(i) && opengl_triangulated_surface(i)){
            glLoadName(Value(i));opengl_triangulated_surface(i)->Display();}
        glPopName();
        glPopName();}
    if(draw_tetrahedralized_volume){
        glPushName(1);
        glPushName(0);
        for(int i=0;i<opengl_tetrahedralized_volume.Size();i++) if(draw_object(i) && opengl_tetrahedralized_volume(i)){
            glLoadName(Value(i));opengl_tetrahedralized_volume(i)->Display();}
        glPopName();
        glPopName();}
    if(draw_implicit_surface){
        glPushName(2);glPushName(0);
        int levelset_count=0;
        for(int i=0;i<opengl_levelset.Size();i++) if(draw_object(i)){
            glLoadName(Value(i));
            if(opengl_levelset(i)){
                if(++levelset_count>50){
                    PHYSBAM_WARNING("Refusing to draw more than 10 levelsets to save memory.");
                    OPENGL_WORLD<T>::Singleton()->Add_String("WARNING: Refusing to draw more than 10 levelsets to save memory.");
                    break;}
                opengl_levelset(i)->Update();
                opengl_levelset(i)->Display();}}
        glPopName();glPopName();}
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
            if(mode==GL_SELECT){
                glPushName(3);
                glPushName(0);
                glPushAttrib(GL_POINT_BIT);
                glPointSize(OPENGL_PREFERENCES::selection_point_size);}
            for(int i=0;i<articulation_points.m;i++){
                glLoadName(i);
                OPENGL_SHAPES::Draw_Dot(articulation_points(i),articulation_point_color,5);}
            OPENGL_SHAPES::Draw_Dot(projected_COM,com_color,10);
            if(mode!=GL_SELECT) for(int i=0;i<articulation_points.m;i+=2){
                OPENGL_SHAPES::Draw_Segment(articulation_points(i),articulation_points(i+1),segment_color,5);}
            if(mode==GL_SELECT){
                glPopName();
                glPopName();
                glPopAttrib();}
            if(mode!=GL_SELECT && selected_joint_id>=0){
                OPENGL_SELECTION::Draw_Highlighted_Vertex(articulation_points(selected_joint_id));}}}

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
            std::string label=LOG::sprintf("F=%.3f %.3f %.3f, T=%.3f %.3f %.3f",forces_and_torques(i).x.x,forces_and_torques(i).x.y,forces_and_torques(i).x.z,forces_and_torques(i).y.x,forces_and_torques(i).y.y,forces_and_torques(i).y.z);
            OpenGL_String(rigid_body_collection.rigid_body_particles.frame(i).t+scale*forces_and_torques(i).x,label);}
        glPopAttrib();}

    if(slice && slice->Is_Slice_Mode()) glPopAttrib();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    int selected_body=selected_surface>=0?selected_surface:selected_volume;
    if(selected_body>=0){
        output_stream<<"Rigid body "<<selected_body<<std::endl;

        // handle case of body no longer being active
        if(!rigid_body_collection.Is_Active(selected_body) || !opengl_triangulated_surface(selected_body)){
            output_stream<<"INACTIVE"<<std::endl;
            return;}

        RIGID_BODY<TV> *body=const_cast<RIGID_BODY<TV>*>(&rigid_body_collection.Rigid_Body(selected_body));
        body->Update_Angular_Velocity();

        if(!body->name.empty()){output_stream<<"Name ="<<body->name<<std::endl;}
        output_stream<<"Mass = "<<body->Mass()<<std::endl;
        output_stream<<"Inertia tensor = "<<body->Inertia_Tensor()<<std::endl;
        output_stream<<"Angular momentum = "<<body->Angular_Momentum()<<std::endl;
        output_stream<<"Kinetic energy = "<<body->Kinetic_Energy()<<std::endl;
        output_stream<<std::endl;

        rigid_body_collection.rigid_body_particles.Print(output_stream,selected_body);

        MATRIX<T,4> body_transform=body->Frame().Matrix();

        if(selected_surface>=0){
            opengl_triangulated_surface(selected_surface)->Print_Selection_Info(output_stream,&body_transform);
            int vertex=opengl_triangulated_surface(selected_surface)->selected_vertex;
            if(vertex>=0){
                TV position=body->simplicial_object->particles.X(vertex);
                output_stream<<"Pointwise velocity = "<<body->Pointwise_Object_Velocity(body->World_Space_Point(position))<<std::endl;}}
        if(selected_volume>=0)
            opengl_tetrahedralized_volume(selected_volume)->Print_Selection_Info(output_stream,&body_transform);}

    else if(selected_joint_id>=0){
        int articulation_id=selected_joint_id;
        int joint_index=(articulation_id+1/2);
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
namespace PhysBAM{
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<float>;
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<double>;
}

