//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <Rigids/Joints/JOINT.h>
#include <Rigids/Joints/JOINT_MESH.h>
#include <Rigids/Muscles/ATTACHMENT_POINT.h>
#include <Rigids/Muscles/MUSCLE.h>
#include <Rigids/Muscles/MUSCLE_LIST.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <OpenGL/OpenGL/OPENGL_AXES.h>
#include <OpenGL/OpenGL/OPENGL_LEVELSET_2D.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D(STREAM_TYPE stream_type,const std::string& basedir_input)
    :OPENGL_COMPONENT<T>(stream_type,"Rigid Geometry Collection 2D"),basedir(basedir_input),frame_loaded(-1),valid(false),show_object_names(false),output_positions(true),draw_velocity_vectors(false),
    draw_individual_axes(false),draw_node_velocity_vectors(false),draw_segmented_curve(true),draw_triangulated_area(false),draw_implicit_curve(false),
    draw_articulation_points(false),draw_forces_and_torques(false),draw_linear_muscles(false),need_destroy_rigid_body_collection(true),selected_curve(-1),selected_area(-1),selected_joint_id(-1),selected_muscle_id(-1),
rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(0)),articulated_rigid_body(0),
velocity_field(stream_type,velocity_vectors,positions,OPENGL_COLOR::Cyan(),0.25,true,true),node_velocity_field(stream_type,node_velocity_vectors,node_positions,OPENGL_COLOR::Magenta(),0.25,true,true)
{
    viewer_callbacks.Set("toggle_velocity_vectors",{[this](){Toggle_Velocity_Vectors();},"Toggle velocity vectors"});
    viewer_callbacks.Set("toggle_individual_axes",{[this](){Toggle_Individual_Axes();},"Toggle individual axes"});
    viewer_callbacks.Set("toggle_output_positions",{[this](){Toggle_Output_Positions();},"Toggle output positions"});
    viewer_callbacks.Set("toggle_show_object_names",{[this](){Toggle_Show_Object_Names();},"Toggle show object names"});
    viewer_callbacks.Set("toggle_node_velocity_vectors",{[this](){Toggle_Node_Velocity_Vectors();},"Toggle node velocity vectors"});
    viewer_callbacks.Set("toggle_draw_mode",{[this](){Toggle_Draw_Mode();},"Toggle draw mode"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_articulation_points",{[this](){Toggle_Articulation_Points();},"Toggle articulation points"});
    viewer_callbacks.Set("toggle_linear_muscles",{[this](){Toggle_Linear_Muscles();},"Toggle linear muscles"});
    viewer_callbacks.Set("toggle_forces_and_torques",{[this](){Toggle_Forces_And_Torques();},"Toggle forces and torques"});

    is_animation=true;
    has_init_destroy_information=true;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D(STREAM_TYPE stream_type,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const std::string& basedir_input)
    :OPENGL_COMPONENT<T>(stream_type,"Rigid Geometry Collection 2D"),basedir(basedir_input),frame_loaded(-1),valid(false),show_object_names(false),output_positions(true),draw_velocity_vectors(false),
    draw_individual_axes(false),draw_node_velocity_vectors(false),draw_segmented_curve(true),draw_triangulated_area(false),draw_implicit_curve(false),
    draw_articulation_points(false),draw_forces_and_torques(false),draw_linear_muscles(false),need_destroy_rigid_body_collection(false),selected_curve(-1),selected_area(-1),selected_joint_id(-1),selected_muscle_id(-1),
    rigid_body_collection(rigid_body_collection),articulated_rigid_body(0),
    velocity_field(stream_type,velocity_vectors,positions,OPENGL_COLOR::Cyan(),0.25,true,true),node_velocity_field(stream_type,node_velocity_vectors,node_positions,OPENGL_COLOR::Magenta(),0.25,true,true)
{
    viewer_callbacks.Set("toggle_velocity_vectors",{[this](){Toggle_Velocity_Vectors();},"Toggle velocity vectors"});
    viewer_callbacks.Set("toggle_individual_axes",{[this](){Toggle_Individual_Axes();},"Toggle individual axes"});
    viewer_callbacks.Set("toggle_output_positions",{[this](){Toggle_Output_Positions();},"Toggle output positions"});
    viewer_callbacks.Set("toggle_show_object_names",{[this](){Toggle_Show_Object_Names();},"Toggle show object names"});
    viewer_callbacks.Set("toggle_node_velocity_vectors",{[this](){Toggle_Node_Velocity_Vectors();},"Toggle node velocity vectors"});
    viewer_callbacks.Set("toggle_draw_mode",{[this](){Toggle_Draw_Mode();},"Toggle draw mode"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_articulation_points",{[this](){Toggle_Articulation_Points();},"Toggle articulation points"});
    viewer_callbacks.Set("toggle_linear_muscles",{[this](){Toggle_Linear_Muscles();},"Toggle linear muscles"});
    viewer_callbacks.Set("toggle_forces_and_torques",{[this](){Toggle_Forces_And_Torques();},"Toggle forces and torques"});

    is_animation=true;
    has_init_destroy_information=true;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
~OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D()
{
    if(need_destroy_rigid_body_collection) delete &rigid_body_collection;
}
//#####################################################################
// Function Set_Draw_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Set_Draw_Mode(const int mode)
{
    draw_segmented_curve=(mode&0x01)!=0;
    draw_triangulated_area=(mode&0x02)!=0;
    draw_implicit_curve=(mode&0x04)!=0;
}
//#####################################################################
// Function Get_Draw_Mode
//#####################################################################
template<class T> int OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Get_Draw_Mode() const
{
    return 0x04*(int)(draw_implicit_curve)+0x02*(int)(draw_triangulated_area)+0x01*(int)(draw_segmented_curve);
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Reinitialize(const bool force,const bool read_geometry)
{
    if(draw && (force || (is_animation && (frame_loaded!=frame)) || (!is_animation && (frame_loaded<0)))){
        valid=false;
        if(!FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/rigid_body_particles",basedir.c_str(),frame))) return;

        rigid_body_collection.Read(stream_type,basedir,frame,&needs_init,&needs_destroy); // TODO: avoiding reading triangulated areas
        if(has_init_destroy_information) for(int i=0;i<needs_destroy.m;i++) Destroy_Geometry(needs_destroy(i));

        std::string arb_state_file=LOG::sprintf("%s/%d/arb_state",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(arb_state_file)){
            if(!articulated_rigid_body) articulated_rigid_body=new ARTICULATED_RIGID_BODY<TV>(rigid_body_collection); // TODO: read in the actual particles
            articulated_rigid_body->Read(stream_type,basedir,frame);}
        else{delete articulated_rigid_body;articulated_rigid_body=0;}

        if(FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/arb_info",basedir.c_str(),frame)))
            Read_Articulated_Information(LOG::sprintf("%s/%d/arb_info",basedir.c_str(),frame));

        std::string filename=LOG::sprintf("%s/%d/rigid_body_forces_and_torques",basedir.c_str(),frame);
        if(FILE_UTILITIES::File_Exists(filename)) FILE_UTILITIES::Read_From_File(stream_type,filename,forces_and_torques);
        else forces_and_torques.Resize(0);

        int max_number_of_bodies(max(opengl_segmented_curve.Size(),rigid_body_collection.rigid_body_particles.Size()));
        // only enlarge array as we read in more geometry to memory
        opengl_segmented_curve.Resize(max_number_of_bodies);
        opengl_triangulated_area.Resize(max_number_of_bodies);
        opengl_levelset.Resize(max_number_of_bodies);
        extra_components.Resize(max_number_of_bodies);
        opengl_axes.Resize(max_number_of_bodies);
        draw_object.Resize(max_number_of_bodies);draw_object.Fill(true);
        use_object_bounding_box.Resize(max_number_of_bodies);use_object_bounding_box.Fill(true);

        // Initialize bodies which have become active
        if(has_init_destroy_information) for(int i=0;i<needs_init.m;i++){
            int id=needs_init(i);PHYSBAM_ASSERT(rigid_body_collection.Is_Active(id));
            Create_Geometry(id);}

        // Only display real bodies (not ghost bodies)
        if(FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/partition",basedir.c_str(),frame))) {
            ARRAY<int> particles_of_this_partition;
            FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/partition",basedir.c_str(),frame),particles_of_this_partition);
            for(int i=0;i<max_number_of_bodies;i++)
                draw_object(i)=false;
            for(int i=0;i<particles_of_this_partition.Size();i++)
                draw_object(particles_of_this_partition(i))=true;}

        // Update active bodies / remove inactive bodies
        for(int id=0;id<rigid_body_collection.rigid_body_particles.Size();id++){
            if(rigid_body_collection.Is_Active(id)) Update_Geometry(id);
            else Destroy_Geometry(id);}
        for(int id=rigid_body_collection.rigid_body_particles.Size();id<opengl_segmented_curve.Size();id++) Destroy_Geometry(id);
        if(FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/colors",basedir.c_str(),frame)))
            FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/colors",basedir.c_str(),frame),colors);
        for(int id=0;id<colors.m;id++){
            if(colors(id)==0) Set_Object_Color(id,OPENGL_COLOR::Green());
            if(colors(id)==1) Set_Object_Color(id,OPENGL_COLOR::Magenta());}

        frame_loaded=frame;
        valid=true;}

    Update_Object_Labels();
}
//#####################################################################
// Function Create_Geometry
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Create_Geometry(const int id)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(id);

    if(!opengl_axes(id)){opengl_axes(id)=new OPENGL_AXES<T>(stream_type);}
    if(rigid_body.simplicial_object && !opengl_segmented_curve(id)){
        opengl_segmented_curve(id)=new OPENGL_SEGMENTED_CURVE_2D<T>(stream_type,*rigid_body.simplicial_object);
        opengl_segmented_curve(id)->draw_vertices=true;
        opengl_segmented_curve(id)->Enslave_Transform_To(*opengl_axes(id));}

    // add triangulated area
    TRIANGULATED_AREA<T>* triangulated_area=rigid_body.template Find_Structure<TRIANGULATED_AREA<T>*>();
    if(triangulated_area && !opengl_triangulated_area(id)){
        triangulated_area->mesh.Initialize_Segment_Mesh();
        opengl_triangulated_area(id)=new OPENGL_TRIANGULATED_AREA<T>(stream_type,*triangulated_area);
        opengl_triangulated_area(id)->Enslave_Transform_To(*opengl_axes(id));}

    if(rigid_body.implicit_object && !opengl_levelset(id)){ // ASSUMES LEVELSET_IMPLICIT_CURVE!
        IMPLICIT_OBJECT<TV>* object_space_implicit_object=rigid_body.implicit_object->object_space_implicit_object;
        if(LEVELSET_IMPLICIT_OBJECT<TV>* levelset_implicit_object=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(object_space_implicit_object)){
            opengl_levelset(id)=new OPENGL_LEVELSET_2D<T>(stream_type,levelset_implicit_object->levelset);
            // TODO: fix me
            //opengl_levelset(id)->Set_Color_Mode(OPENGL_LEVELSET_2D<T>::COLOR_GRADIENT);
            opengl_levelset(id)->draw_cells=true;opengl_levelset(id)->draw_curve=false;opengl_levelset(id)->draw_area=false;
            opengl_levelset(id)->Enslave_Transform_To(*opengl_axes(id));}}

    // add extra components
    if(opengl_triangulated_area(id)){
        std::string filename_pattern=LOG::sprintf("%s/accumulated_impulses_%d.%%d",basedir.c_str(),id);
        if(FILE_UTILITIES::Frame_File_Exists(filename_pattern,frame)){
            LOG::cout<<"Adding accumulated impulses to rigid body "<<id<<std::endl;
            OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>* component=
                new OPENGL_COMPONENT_TRIANGULATED_AREA_BASED_VECTOR_FIELD<T>(stream_type,*triangulated_area,filename_pattern);
            component->opengl_vector_field.Enslave_Transform_To(*opengl_axes(id));
            extra_components(id).Append(component);}}
}
//#####################################################################
// Function Update_Geometry
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Update_Geometry(const int id)
{
    if(opengl_axes(id)) *opengl_axes(id)->frame=Convert_2d_To_3d(rigid_body_collection.Rigid_Body(id).Frame());
    if(opengl_triangulated_area(id)){
        std::string color_map_filename=LOG::sprintf("%s/%d/stress_map_of_triangulated_area_%d",basedir.c_str(),frame,id);
        if(FILE_UTILITIES::File_Exists(color_map_filename)){
            if(!opengl_triangulated_area(id)->color_map) opengl_triangulated_area(id)->color_map=new ARRAY<OPENGL_COLOR>;
            FILE_UTILITIES::Read_From_File(stream_type,color_map_filename,*opengl_triangulated_area(id)->color_map);}
        else if(opengl_triangulated_area(id)->color_map){delete opengl_triangulated_area(id)->color_map;opengl_triangulated_area(id)->color_map=0;}}
    for(int i=0;i<extra_components(id).m;i++) extra_components(id)(i)->Set_Frame(frame);
}
//#####################################################################
// Function Destroy_Geometry
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Destroy_Geometry(const int id)
{
    if(opengl_segmented_curve(id)){delete opengl_segmented_curve(id);opengl_segmented_curve(id)=0;}
    if(opengl_triangulated_area(id)){delete opengl_triangulated_area(id)->color_map;delete opengl_triangulated_area(id);opengl_triangulated_area(id)=0;}
    if(opengl_levelset(id)){delete opengl_levelset(id);opengl_levelset(id)=0;}
    if(opengl_axes(id)){delete opengl_axes(id);opengl_axes(id)=0;}
    extra_components(id).Delete_Pointers_And_Clean_Memory();
    draw_object(id)=false;
}
//#####################################################################
// Function Update_Object_Labels
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Update_Object_Labels()
{
    int number_of_drawn_bodies=draw_object.Count_Matches(true);
    if(draw_velocity_vectors){
        positions.Resize(number_of_drawn_bodies);
        velocity_vectors.Resize(number_of_drawn_bodies);}
    if(draw_node_velocity_vectors){
        node_positions.Resize(number_of_drawn_bodies*4);
        node_velocity_vectors.Resize(number_of_drawn_bodies*4);}

    int idx=0;
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
        if(draw_object(i)){
            if(draw_velocity_vectors || draw_node_velocity_vectors){idx++;
                if(draw_velocity_vectors){
                    positions(idx)=rigid_body_collection.rigid_body_particles.frame(i).t;
                    velocity_vectors(idx)=rigid_body_collection.rigid_body_particles.twist(i).linear;}
                if(draw_node_velocity_vectors){
                    //only valid for squares . . .
                    RIGID_BODY<TV>* rigid_body=&rigid_body_collection.Rigid_Body(i);

                    node_positions((idx-1)*4+1)=rigid_body->World_Space_Point(VECTOR<T,2>(1,1));
                    node_velocity_vectors((idx-1)*4+1)=rigid_body->Pointwise_Object_Velocity(rigid_body->World_Space_Vector(VECTOR<T,2>(1,1)));

                    node_positions((idx-1)*4+2)=rigid_body->World_Space_Point(VECTOR<T,2>(1,-1));
                    node_velocity_vectors((idx-1)*4+2)=rigid_body->Pointwise_Object_Velocity(rigid_body->World_Space_Vector(VECTOR<T,2>(1,-1)));

                    node_positions((idx-1)*4+3)=rigid_body->World_Space_Point(VECTOR<T,2>(-1,1));
                    node_velocity_vectors((idx-1)*4+3)=rigid_body->Pointwise_Object_Velocity(rigid_body->World_Space_Vector(VECTOR<T,2>(-1,1)));

                    node_positions((idx-1)*4+4)=rigid_body->World_Space_Point(VECTOR<T,2>(-1,-1));
                    node_velocity_vectors((idx-1)*4+4)=rigid_body->Pointwise_Object_Velocity(rigid_body->World_Space_Vector(VECTOR<T,2>(-1,-1)));}}
            if(opengl_segmented_curve(i)){
                if(output_positions){
                    opengl_segmented_curve(i)->Set_Name(LOG::sprintf("%s <%.3f %.3f> [w=%.3f]",rigid_body_collection.Rigid_Body(i).name.c_str(),rigid_body_collection.rigid_body_particles.frame(i).t.x,rigid_body_collection.rigid_body_particles.frame(i).t.y,rigid_body_collection.rigid_body_particles.twist(i).angular.x));}
                else opengl_segmented_curve(i)->Set_Name(rigid_body_collection.Rigid_Body(i).name);}}}
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
        if(draw_object(i)){
            if(opengl_segmented_curve(i)){
                if(output_positions){
                rigid_body_collection.Rigid_Body(i).Update_Angular_Velocity();}}}}
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::File_Exists(LOG::sprintf("%s/%d/rigid_body_particles",basedir.c_str(),frame_input));
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
}
//#####################################################################
// Function Use_Bounding_Box
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Use_Bounding_Box() const
{
    int num_drawable_objects=0;
    for(int i=0;i<opengl_segmented_curve.Size();i++)
        if(draw_object(i) && use_object_bounding_box(i))
            num_drawable_objects++;
    return (draw && num_drawable_objects>0);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Bounding_Box() const
{
    RANGE<VECTOR<T,3> > box;
    if(draw)
        for(int i=0;i<opengl_segmented_curve.Size();i++)
            if(draw_object(i) && use_object_bounding_box(i) && opengl_segmented_curve(i)){
                box=RANGE<VECTOR<T,3> >::Combine(box,opengl_segmented_curve(i)->Bounding_Box());}
    return box;
}
//#####################################################################
// Function Get_Selection_Priority
//#####################################################################
template<class T> int OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Get_Selection_Priority(ARRAY_VIEW<GLuint> indices)
{
    if(!indices.m) return -1;
    const static int priority[]={0,0,95,95};
    if(indices(0)==0) // segmented curve
        return opengl_segmented_curve(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
    if(indices(0)==1) // triangulated area
        return opengl_triangulated_area(indices(1))->Get_Selection_Priority(indices.Array_View(2,indices.m-2));
    return priority[indices(0)];
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers)
{
    if(indices(0)==0 && (modifiers&GLUT_ACTIVE_CTRL)){
        draw_object(indices(1))=false;
        return false;}
    if(indices(0)==0){
        selected_curve=indices(1);
        return opengl_segmented_curve(indices(1))->Set_Selection(indices.Array_View(2,indices.m-2),modifiers);}
    else if(indices(0)==1){
        selected_area=indices(1);
        return opengl_triangulated_area(indices(1))->Set_Selection(indices.Array_View(2,indices.m-2),modifiers);}
    else if(indices(0)==4) // articulation joints
        selected_joint_id=indices(1);
    else if(indices(0)==5)
        selected_muscle_id=indices(1);
    else return false;
    return true;
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Clear_Selection()
{
    if(selected_curve>=0){
        opengl_segmented_curve(selected_curve)->Clear_Selection();
        selected_curve=-1;}
    if(selected_area>=0){
        opengl_triangulated_area(selected_area)->Clear_Selection();
        selected_area=-1;}
    selected_joint_id=-1;
    selected_muscle_id=-1;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Print_Selection_Info(std::ostream &output_stream) const
{
    int selected_body=selected_curve>=0?selected_curve:selected_area;
    if(selected_body>=0){
        output_stream<<"Rigid body "<<selected_body<<std::endl;

        // handle case of body no longer being active
        if(!rigid_body_collection.Is_Active(selected_body) || !opengl_segmented_curve(selected_body)){
            output_stream<<"INACTIVE"<<std::endl;
            return;}

        RIGID_BODY<TV> *body=const_cast<RIGID_BODY<TV>*>(&rigid_body_collection.Rigid_Body(selected_body));
        //body->Update_Angular_Velocity();

        if(!body->name.empty()){output_stream<<"Name = "<<body->name<<std::endl;}
        rigid_body_collection.rigid_body_particles.Print(output_stream,selected_body);

        MATRIX<T,3> body_transform=body->Frame().Matrix();

        if(selected_curve>=0){
            opengl_segmented_curve(selected_curve)->Print_Selection_Info(output_stream,&body_transform);
            int vertex=opengl_segmented_curve(selected_curve)->selected_vertex;
            if(vertex>=0){
                TV position=body->simplicial_object->particles.X(vertex);
                output_stream<<"Pointwise velocity = "<<body->Pointwise_Object_Velocity(body->World_Space_Point(position))<<std::endl;}}
        if(selected_area>=0)
            opengl_triangulated_area(selected_area)->Print_Selection_Info(output_stream,&body_transform);}

    else if(selected_joint_id>=0){
        JOINT_ID joint_id(selected_joint_id/2);
        JOINT<TV>* joint=articulated_rigid_body->joint_mesh(joint_id);
        const RIGID_BODY<TV>* parent=articulated_rigid_body->Parent(joint_id),*child=articulated_rigid_body->Child(joint_id);
        FRAME<TV> current_frame=joint->Compute_Current_Joint_Frame(*parent,*child);
        output_stream<<"Joint id "<<joint_id<<" ("<<(!joint->name.empty()?joint->name:"UNNAMED")<<")"<<std::endl;
        output_stream<<"Frame = "<<(selected_joint_id%2==0?"child":"parent")<<std::endl;
        //output_stream<<"Articulation point = "<<articulation_points(selected_joint_id)<<std::endl;
        VECTOR<T,2> ap1=parent->World_Space_Point(joint->F_pj().t),ap2=child->World_Space_Point(joint->F_cj().t),location=(T).5*(ap1+ap2);
        output_stream<<"Parent = "<<parent->name<<std::endl;
        output_stream<<"Child = "<<child->name<<std::endl;
        output_stream<<"Joint translation = "<<current_frame.t<<std::endl;
        output_stream<<"Joint rotation angle = "<<current_frame.r.Angle()<<std::endl;

        VECTOR<T,2> current_relative_velocity=-RIGID_BODY<TV>::Relative_Velocity(*parent,*child,location); // child w.r.t. parent!
        VECTOR<T,1> current_relative_angular_velocity=-RIGID_BODY<TV>::Relative_Angular_Velocity(*parent,*child); // child w.r.t. parent!
        output_stream<<"Relative velocity at joint = "<<current_relative_velocity<<std::endl;
        output_stream<<"Relative angular velocity = "<<current_relative_angular_velocity<<std::endl;}
    else if(selected_muscle_id>=0){
        int muscle_id=selected_muscle_id;
        MUSCLE<TV>& muscle=*articulated_rigid_body->muscle_list->muscles(muscle_id);
        ATTACHMENT_POINT<TV>* attachment_point_1=muscle.attachment_point_1;
        ATTACHMENT_POINT<TV>* attachment_point_2=muscle.attachment_point_2;
        output_stream<<"Muscle "<<muscle_id<<" ("<<(!muscle.name.empty()?muscle.name:"UNNAMED")<<")"<<std::endl;
        output_stream<<"Optimal length = "<<muscle.optimal_length<<std::endl;
        output_stream<<"Peak force = "<<muscle.peak_force<<std::endl;
        output_stream<<"Pennation angle = "<<muscle.pennation_angle<<std::endl;
        output_stream<<"Tendon slack length = "<<muscle.tendon_slack_length<<std::endl;
        output_stream<<"Maximum shortening velocity = "<<muscle.max_shortening_velocity<<std::endl;
        output_stream<<"Origin = ("<<attachment_point_1->rigid_body.name<<", "<<attachment_point_1->object_space_position<<")"<<std::endl;
        output_stream<<"Insertion = ("<<attachment_point_2->rigid_body.name<<", "<<attachment_point_2->object_space_position<<")"<<std::endl;
        if(muscle.via_points.m){
            output_stream<<"Via points: ";
            for(int i=0;i<muscle.via_points.m;i++)
                output_stream<<"("<<muscle.via_points(i)->rigid_body.name<<", "<<muscle.via_points(i)->object_space_position<<") ";
            output_stream<<std::endl;}
        output_stream<<std::endl;
        output_stream<<"Total length = "<<muscle.Total_Length()<<std::endl;
        output_stream<<"Total velocity = "<<muscle.Total_Velocity()<<std::endl;
        if(articulated_rigid_body->muscle_activations.Valid_Index(muscle_id)) output_stream<<"Activation = "<<articulated_rigid_body->muscle_activations(muscle_id)<<std::endl;}
}
//#####################################################################
// Function Set_Draw_Object
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Set_Draw_Object(int i,bool draw_it)
{
    if(i>draw_object.Size()) draw_object.Resize(i);
    draw_object(i)=draw_it;
}
//#####################################################################
// Function Set_Object_Color
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Set_Object_Color(int i,const OPENGL_COLOR &color_input)
{
    if(!opengl_segmented_curve(i)) return;
    opengl_segmented_curve(i)->color=color_input;
    opengl_segmented_curve(i)->vertex_color=color_input;
}
//#####################################################################
// Function Get_Draw_Object
//#####################################################################
template<class T> bool OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Get_Draw_Object(int i) const
{
    return draw_object(i);
}
//#####################################################################
// Function Set_Use_Object_Bounding_Box
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Set_Use_Object_Bounding_Box(int i,bool use_it)
{
    use_object_bounding_box(i)=use_it;
}
//#####################################################################
// Function Set_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Set_Vector_Size(T size)
{
    velocity_field.size=size;
    node_velocity_field.size=size;
}
//#####################################################################
// Function Toggle_Velocity_Vectors
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Toggle_Velocity_Vectors()
{
    draw_velocity_vectors=!draw_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Individual_Axes
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Toggle_Individual_Axes()
{
    draw_individual_axes=!draw_individual_axes;
    Reinitialize();
}
//#####################################################################
// Function Toggle_Output_Positions
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Toggle_Output_Positions()
{
    output_positions=!output_positions;
    Update_Object_Labels();
}
//#####################################################################
// Function Toggle_Show_Object_Names
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Toggle_Show_Object_Names()
{
    show_object_names=!show_object_names;
}
//#####################################################################
// Function Toggle_Node_Velocity_Vectors
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Toggle_Node_Velocity_Vectors()
{
    draw_node_velocity_vectors=!draw_node_velocity_vectors;
    Reinitialize();
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Increase_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1.1);
    node_velocity_field.Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Decrease_Vector_Size()
{
    velocity_field.Scale_Vector_Size(1/1.1);
    node_velocity_field.Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Draw_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Toggle_Draw_Mode()
{
    int mode=Get_Draw_Mode();
    mode=(mode+1)%8;
    Set_Draw_Mode(mode);
}
//#####################################################################
// Selection Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Selection_Bounding_Box() const
{
    if(selected_curve>=0) return opengl_segmented_curve(selected_curve)->Selection_Bounding_Box();
    if(selected_area>=0) return opengl_triangulated_area(selected_area)->Selection_Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Read_Articulated_Information
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Read_Articulated_Information(const std::string& filename)
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    TYPED_ISTREAM typed_input(*input,stream_type);
    // this will need to be changed to reflect multiple articulation points per rigid body
    int numpoints=0;Read_Binary(typed_input,numpoints);articulation_points.Exact_Resize(numpoints);
    for(int i=0;i<numpoints;i++) Read_Binary(typed_input,articulation_points(i));
    delete input;
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Display() const
{
    if(draw){
        glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);

        GLint mode=0;
        glGetIntegerv(GL_RENDER_MODE,&mode);

        if(draw_segmented_curve){
            glPushName(0);
            glPushName(0);
            for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
                glLoadName(Value(i));
                if(draw_object(i) && opengl_segmented_curve(i)) opengl_segmented_curve(i)->Display();}
            glPopName();
            glPopName();}
        if(draw_triangulated_area){
            glPushName(1);
            glPushName(0);
            for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
                glLoadName(Value(i));
                if(draw_object(i) && opengl_triangulated_area(i)) opengl_triangulated_area(i)->Display();}
            glPopName();
            glPopName();}
        if(draw_implicit_curve){
            glPushName(2);
            glPushName(0);
            for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++){
                glLoadName(Value(i));
                if(draw_object(i) && opengl_levelset(i)) opengl_levelset(i)->Display();}
            glPopName();
            glPopName();}
        if(draw_individual_axes)
            for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++)
                if(draw_object(i) && opengl_axes(i)) opengl_axes(i)->Display();

        // Articulated rigid bodies
        if(articulated_rigid_body){
            if(draw_articulation_points){
                OPENGL_COLOR articulation_point_color(0.5,0.5,0.5),segment_color(0,0,1);
                if(mode==GL_SELECT){
                    glPushName(3);
                    glPushName(0);
                    glPushAttrib(GL_POINT_BIT);
                    glPointSize(OPENGL_PREFERENCES::selection_point_size);}
                for(int i=0;i<articulation_points.m;i++){
                    glLoadName(i);
                    OPENGL_SHAPES::Draw_Dot(articulation_points(i),articulation_point_color,5);}
                if(mode!=GL_SELECT)
                    for(int i=0;i<articulation_points.m;i+=2)
                        OPENGL_SHAPES::Draw_Segment(articulation_points(i),articulation_points(i+1),segment_color,5);
                if(mode==GL_SELECT){
                    glPopName();
                    glPopName();
                    glPopAttrib();}
                if(mode!=GL_SELECT && selected_joint_id)
                    OPENGL_SELECTION::Draw_Highlighted_Vertex(articulation_points(selected_joint_id));}

            if(draw_linear_muscles){
                glPushAttrib(GL_LINE_BIT|GL_ENABLE_BIT|GL_CURRENT_BIT);
                // draw all muscles
                OPENGL_COLOR muscle_color(1,0,0);muscle_color.Send_To_GL_Pipeline();
                glLineWidth(mode==GL_SELECT?OPENGL_PREFERENCES::selection_line_width:5);
                glPushName(4);
                glPushName(0);
                for(int i=0;i<articulated_rigid_body->muscle_list->muscles.m;i++){
                    MUSCLE<TV>& muscle=*articulated_rigid_body->muscle_list->muscles(i);
                    glLoadName(i);OpenGL_Begin(GL_LINE_STRIP);
                    OpenGL_Vertex(muscle.attachment_point_1->Embedded_Position());
                    for(int t=0;t<muscle.via_points.m;t++) OpenGL_Vertex(muscle.via_points(t)->Embedded_Position());
                    OpenGL_Vertex(muscle.attachment_point_2->Embedded_Position());
                    OpenGL_End();}
                glPopName();
                glPopName();
                // highligh selected one
                if(mode!=GL_SELECT && selected_muscle_id){
                    MUSCLE<TV>& muscle=*articulated_rigid_body->muscle_list->muscles(selected_muscle_id);
                    glLineWidth(OPENGL_PREFERENCES::highlighted_line_width);OPENGL_PREFERENCES::selection_highlight_color.Send_To_GL_Pipeline();
                    OpenGL_Begin(GL_LINE_STRIP);
                    OpenGL_Vertex(muscle.attachment_point_1->Embedded_Position());
                    for(int t=0;t<muscle.via_points.m;t++) OpenGL_Vertex(muscle.via_points(t)->Embedded_Position());
                    OpenGL_Vertex(muscle.attachment_point_2->Embedded_Position());
                    OpenGL_End();}
                glPopAttrib();}}


        if(mode!=GL_SELECT)
        {
            if(draw_velocity_vectors) velocity_field.Display();
            if(draw_node_velocity_vectors) node_velocity_field.Display();

            if(draw_forces_and_torques && forces_and_torques.Size()==rigid_body_collection.rigid_body_particles.Size()){
                T scale=(T)velocity_field.size/24;
                OPENGL_COLOR::Yellow().Send_To_GL_Pipeline();
                OpenGL_Begin(GL_LINES);
                for(int i=0;i<forces_and_torques.Size();i++)
                    OPENGL_SHAPES::Draw_Arrow(rigid_body_collection.rigid_body_particles.frame(i).t,rigid_body_collection.rigid_body_particles.frame(i).t+scale*forces_and_torques(i).x);
                OpenGL_End();
                for(int i=0;i<forces_and_torques.Size();i++){
                    std::string label=LOG::sprintf("F=%.3f %.3f, T=%.3f",forces_and_torques(i).x.x,forces_and_torques(i).x.y,forces_and_torques(i).y);
                    OpenGL_String(rigid_body_collection.rigid_body_particles.frame(i).t+scale*forces_and_torques(i).x,label);}}

            for(int i=0;i<extra_components.Size();i++)
                for(int j=0;j<extra_components(i).m;j++)
                    extra_components(i)(j)->Display();

            if(show_object_names){
                glColor3f(1,1,1);
                for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++)
                    if(draw_object(i) && rigid_body_collection.Rigid_Body(i).name.length())
                        OpenGL_String(rigid_body_collection.rigid_body_particles.frame(i).t,rigid_body_collection.Rigid_Body(i).name);}
        }
        glPopAttrib();}
}
//#####################################################################
// Function Toggle_Articulation_Points
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Toggle_Articulation_Points()
{
    draw_articulation_points=!draw_articulation_points;
}
//#####################################################################
// Function Toggle_Articulation_Points
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Toggle_Linear_Muscles()
{
    draw_linear_muscles=!draw_linear_muscles;
}
//#####################################################################
// Function Toggle_Forces_And_Torques
//#####################################################################
template<class T> void OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>::
Toggle_Forces_And_Torques()
{
    draw_forces_and_torques=!draw_forces_and_torques;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<float>;
template class OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<double>;
}
