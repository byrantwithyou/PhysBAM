//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Unnur Gretarsdottir.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHRINK_WRAP_EXAMPLE
//#####################################################################
#ifndef __INNER_LAYER_SIM_EXAMPLE__
#define __INNER_LAYER_SIM_EXAMPLE__


#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include "Inner_Layer_Sim/SKIN_DEFORMABLE_OBJECT_3D.h"
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
#include <Geometry/TRIANGULATED_SURFACE_LIST.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>

namespace PhysBAM{
template<class T,class RW>
class INNER_LAYER_SIM_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:

    std::string input_file;
    ARRAY<bool> ls_inside;
    ARRAY<VECTOR_3D<T> > ls_normal;
    T skin_ls_stiff, skin_ls_damping, skin_las_stiff, skin_las_damping;
    T skin_bending_stiff,skin_bending_damping,bending_cutoff_ratio;
    int num_skin_particles;
    ARRAY<int> closest_muscle_arr;
    ARRAY<int> closest_triangle_arr;
    ARRAY<VECTOR_3D<T> >barycentric_coordinate_arr;
    SKIN_DEFORMABLE_OBJECT_3D<T> *my_skin;
    T threshold_for_skin_to_rigid_bodies_springs;
    ARRAY<TRIANGULATED_SURFACE<T> *> muscle_tris;
    ARRAY<TRIANGULATED_SURFACE<T> *> muscle_tris_old;
    char muscle_directory[256]; char levelset_dir[256];
    bool skip_big_ones;
    bool write_levelsets;
    bool reset_rest_lengths;
    bool set_v_to_zero_after_snap_to_levelset;
    bool free_border_constrained_nodes;
    bool use_bending;
    bool reset_bending;
    LEVELSET_IMPLICIT_SURFACE<T> *levelset;
    T snap_depth;

    INNER_LAYER_SIM_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE)
    {   
        /////////////////////////////////////////
        //Parameters you could/should change
        //////////////////////////////////////////

        //when inner layer snaps to the levelset, all nodes whose depth < snap_depth 
        //get pushed to the levelset (note nodes inside the levelset have a negative 
        //depth and should always get pushed, so snap_depth should be >= 0. 
        snap_depth = (T).001; 
        
        //the sim does "settling" by running time between the frames. Therefore,
        //increasing the frame rate decreases the "settling" time and vice versa.
        frame_rate = 6; 

        //pretty straightforward - stiffness of springs in the inner layer
        skin_ls_stiff = (T).1*2/(1+sqrt((T)2)); skin_ls_damping = 2;
        skin_las_stiff = (T).1*2*4/(1+sqrt((T)2)); skin_las_damping = 4;
        skin_bending_stiff=(T)1e-3;skin_bending_damping=1e-3;bending_cutoff_ratio=(T).5;
         
        //random sim info and stuff you know about
        restart=false;restart_frame=16; last_frame=10*(int)frame_rate;
        output_directory="Inner_Layer_Sim/output";
        input_file="Inner_Layer_Sim/inner_layer.tri.gz";
        sprintf(muscle_directory, "../quasistatics/Muscle/output/");
        sprintf(levelset_dir, "Inner_Layer_Sim/Levelsets");
        solids_parameters.cfl=(T)5.9;  solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        solids_parameters.use_constant_mass = true; solids_parameters.perform_self_collision=false; solids_parameters.collide_with_interior=false;

        //not really relevant to sim quality
        skip_big_ones = false; //will automatically skip some of the muscles/bones in the chest area in order to make sim faster.
        threshold_for_skin_to_rigid_bodies_springs = 100.0;  //in this sim, you should just keep this # high (above 1.0 should be fine)
         write_output_files = false; //we do outputting manually in Postprocess_Substep
        write_levelsets = true; //write muscle/bone leveslets each frame
        reset_rest_lengths = false; // reset spring restlengths after settling
        set_v_to_zero_after_snap_to_levelset = true; //set the node v to zero after it is snapped to the leveslet
        free_border_constrained_nodes = true; // if a node has a neighbor that is not constrained, don't constrain it
        use_bending = true;
        reset_bending = true; //resets the bending rest angle based on the current mesh
      }

    ~INNER_LAYER_SIM_EXAMPLE()
    {}

//#####################################################################
// Function Postprocess_Frame 
//=======================================================================
// This is the meat of the sim. At each frame, we do the projection and 
// then let the sim run to allow the skin to "settle". Note that for 
// some reason this does not get called at the end of frame 0, so all the 
// numbers are adjusted so stuff runs correctly in spite of that... 
//#####################################################################
void Postprocess_Frame(const int frame){
    //Write out sim data. Note that since this happens 3 times per frame, we adjust the frame # written
    Write_Output_Files_Helper(3*(frame-1));

    //reset the restlengths in the skin mesh at each frame (can be changed)
    if(reset_rest_lengths) Reset_Restlengths();
    if(reset_bending && use_bending) Reset_Bending(false);

    //read in levelset for the next frame (don't change)
    char levelset_file[256]; sprintf(levelset_file, "%s/arm_skin_%d.phi", levelset_dir, frame);
    delete levelset; levelset = new LEVELSET_IMPLICIT_SURFACE<T>(*(new GRID<TV>()), *(new ARRAY<T,VECTOR<int,3> >()));
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(levelset_file);levelset->template Read<RW>(*input);delete input;
    levelset->Compute_Cell_Minimum_And_Maximum();

    //Find the closest point on muscle for each particle in the inner layer, read in the muscles at the next frame,
    //and project the inner layer particles by the amount that their closest muscle point moves. (don't change)
    muscle_tris_old.Resize(muscle_tris.m);
    for(int i=0;i<muscle_tris.m;i++) {muscle_tris_old(i) = muscle_tris(i);}
    Find_Spring_Endpoints();
    Read_Muscle_Tris_From_Quasistatic_Simm_Data(frame);    
    Project_Skin_Particles();
    Write_Output_Files_Helper(3*(frame-1)+1);

    //For each inner layer particle which is inside or close to the muscle levelset, push it to the levelset
    //and constrain it's motion to stay in the plane of the levelset. Note - this can be done here or at each
    //timestep (which is set in Skin_Deformable_Object_3D) (can be changed)
    Snap_To_Levelset();
    Write_Output_Files_Helper(3*(frame-1)+2);

    //set velocities to 0 so their motion from the last frame doesn't affect this frame. (don't change)
    Reset_Velocities();

    //Other cleanup stuff (don't change)
    for(int i=0;i<muscle_tris_old.m;i++) {delete (muscle_tris_old(i));}
    std::cout<<"Frame number "<<frame-1<<" completed"<<std::endl;
}

//#######################################################################
// Function Set_External_Velocities and Zero_Out_Enslaved_Velocity_Nodes
// =======================================================================
// Those nodes that have ls_inside set are constrained to the plane of the levelset
//#######################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time) {
    for(int c=0;c<V.m;c++) {
        if(ls_inside(c)) { 
            V(c)-=VECTOR_3D<T>::Dot_Product(V(c),ls_normal(c))*ls_normal(c); 
        }}
}
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time) {    
    for(int c=0;c<V.m;c++) {
        if(ls_inside(c)) { 
            V(c)-=VECTOR_3D<T>::Dot_Product(V(c),ls_normal(c))*ls_normal(c); 
        }}
}

//#####################################################################
// Snap_To_Levelset
// ==================
// Takes each particle in the inner layer and if it is inside the levelset, or
// outside the levelset, but within some threshold and pushes it to the levelset.
// Those particles which are pushed have ls_inside set to true, so they can later
// be constrained to the plane of the levelset.
//#####################################################################
void Snap_To_Levelset(){
    //Push each node to the levelset
    for(int p=0;p<my_skin->particles.array_collection->Size();p++) {
        VECTOR_3D<T> current_location = my_skin->particles.X(p);T depth;
        bool outside = levelset->Lazy_Outside_Extended_Levelset_And_Value(current_location, depth); 
        bool inside = levelset->Lazy_Inside_Extended_Levelset_And_Value(current_location, depth); 
        VECTOR_3D<T> n = levelset->Extended_Normal(current_location);
        //note - Lazy Inside/Outside does not always return correct bool(??) Use sign of depth1 instead...
        ls_inside(p) = false;
        if(depth <= snap_depth) {
            my_skin->particles.X(p) = current_location - depth*n;
            ls_inside(p) = true;
            if(set_v_to_zero_after_snap_to_levelset) my_skin->particles.V(p) = VECTOR_3D<T>(0,0,0);
            ls_normal(p) = levelset->Extended_Normal(my_skin->particles.X(p));}}

    //"Free" the border nodes (can be turned on/off)
    if(free_border_constrained_nodes) Modify_LS_Inside();
}
//#####################################################################
// Modify_LS_Inside
// ========================
// Looks at each particle in the inner layer and checks if ls_inside is 
// true. ls_inside is true if the particle was pushed to the levelset 
// the last time we did Levelset_Snap (these are the nodes that are constrained
// to move in the plane of the levelset). This function will "free" the nodes
// who have neighbors which are unconstrained. This allows nodes which are on
// the border of regions which are off the levelset to slide off the levleset.
//#####################################################################
void Modify_LS_Inside(){
    ARRAY<bool >ls_inside_new; ls_inside_new.Resize(ls_inside.m);
    for(int i=0;i<ls_inside_new.m;i++) {
        ARRAY<int> neighbors = (*(my_skin->triangulated_surface->triangle_mesh.neighbor_nodes))(i);
        ls_inside_new(i) = ls_inside(i); 
        for(int n=0;n<neighbors.m;n++) {
            if(!ls_inside(neighbors(n))) ls_inside_new(i) = false;
        }
    }
    for(int i=0;i<ls_inside_new.m;i++) ls_inside(i) = ls_inside_new(i);
}

////////////////////////////////////////Stuff you probably don't need to mess with////////////////////////////////////////////////

//#####################################################################
// Function Find_Skin_Endpoints
//====================================
// For each particle in the inner layer, finds the closest point on the muscles.
// Starts searching for close ones and increases the proximitly threshold
// incrementally. In this sim, threshold_for_skin_to_rigid_bodies_springs is 
// set very high, so it'll find a closest point on the muscles for every
// particle in the inner layer.
//#####################################################################
void Find_Spring_Endpoints() {    
    T initial_proximity_threshold; 
    T proximity_threshold_squared; 
    T upper_threshold = threshold_for_skin_to_rigid_bodies_springs;
    bool do_over = false;
    for (int i=1;i<=num_skin_particles;i++) {
        if(!do_over){initial_proximity_threshold = (T)0.004; proximity_threshold_squared = initial_proximity_threshold * initial_proximity_threshold;}
        closest_muscle_arr(i) = 0; T distance_sq = 7000000;
        if(initial_proximity_threshold >= upper_threshold) {do_over=false; continue; std::cout<<"!!!Uh Oh - you should increase the upper_threshold!!!"<<std::endl;}
        VECTOR_3D<T> skin_particle = my_skin->particles.X(i);
        int max=muscle_tris.m; 
        for (int j=1;j<=max;j++) {
            TRIANGULATED_SURFACE<T> *muscle_surface = muscle_tris(j);
            ARRAY<int> nearby_triangles;
            int closest_triangle=0;
            VECTOR_3D<T> closest_point_on_triangle, closest_point_on_surface, barycentric_coordinates, closest_barycentric_coordinates;
            T triangle_point_distance_squared, closest_triangle_point_distance_squared;
            muscle_surface->triangle_hierarchy->Intersection_List(skin_particle,nearby_triangles,initial_proximity_threshold);
            for (int k=1;k<=nearby_triangles.m;k++) {
                closest_point_on_triangle=(*muscle_surface->triangle_list)(nearby_triangles(k)).Closest_Point(skin_particle,barycentric_coordinates);
                triangle_point_distance_squared=(skin_particle-closest_point_on_triangle).Magnitude_Squared();
                if (triangle_point_distance_squared>proximity_threshold_squared) continue;
                if (closest_triangle==0||triangle_point_distance_squared<closest_triangle_point_distance_squared) {
                    closest_triangle=nearby_triangles(k); closest_triangle_point_distance_squared=triangle_point_distance_squared;
                    closest_point_on_surface=closest_point_on_triangle; closest_barycentric_coordinates = barycentric_coordinates;} }
            if (closest_triangle==0) continue;  //this muscle is too far away
            if (closest_triangle_point_distance_squared <= distance_sq) { //this muscle is the new closest one
                closest_muscle_arr(i) = j;
                closest_triangle_arr(i) = closest_triangle;
                barycentric_coordinate_arr(i) = closest_barycentric_coordinates;}}
        do_over=false;
        if (closest_muscle_arr(i) == 0) {i--; initial_proximity_threshold *= 2; proximity_threshold_squared = initial_proximity_threshold * initial_proximity_threshold; do_over=true; }
    }
    std::cout<<"Done Finding Spring Endpoints"<<std::endl;
}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() {
    Get_Initial_Data();
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*my_skin->triangulated_surface;
    my_skin->Add_Edge_Springs(triangulated_surface.triangle_mesh,skin_ls_stiff,skin_ls_damping);
    my_skin->Add_Altitude_Springs(triangulated_surface.triangle_mesh,skin_las_stiff,skin_las_damping);

    //bending can be turned on/off to tune sim. In general, it's good, but it can hose the sim because
    //so many nodes are constrained to the levelset and that's fighting against the bending forces...
    if(use_bending){
    solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements(triangulated_surface.triangle_mesh,skin_bending_stiff,skin_bending_damping);
    solids_parameters.deformable_body_parameters.list(1).bending_elements(1)->Set_Area_Cutoff_From_Triangulated_Surface(triangulated_surface,bending_cutoff_ratio);
    Reset_Bending(false);} //if you use bending, it should probably be set to zero
}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    //Skin Mesh
     my_skin = new SKIN_DEFORMABLE_OBJECT_3D<T>(); my_skin->example = this; my_skin->do_regular_collisions = false;
    solids_parameters.deformable_body_parameters.list.deformable_objects.Append(my_skin);
    int index = solids_parameters.deformable_body_parameters.list.deformable_objects.m;
    solids_parameters.deformable_body_parameters.list.deformable_objects(index)->Allocate_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).triangulated_surface;
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_file);triangulated_surface.template Read<RW>(*input);delete input;
    std::cout << "total vertices = " << particles.array_size << std::endl;std::cout << "total triangles = " << triangle_mesh.triangles.m << std::endl;
    particles.store_velocity = false; particles.Update_Velocity();  // forcing velocity resize - not sure why this is needed
    particles.store_mass = false; particles.Store_Mass(); // forcing mass resize - not sure why
    triangulated_surface.Set_Density(1000);triangulated_surface.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    triangulated_surface.Update_Bounding_Box();    triangulated_surface.triangle_mesh.Initialize_Neighbor_Nodes();

    //Resize Arrays
    num_skin_particles = particles.array_collection->Size(); for(int p=0;p<num_skin_particles;p++) {particles.mass(p) = 0.0;}
    closest_muscle_arr.Resize(num_skin_particles); for(int p=0;p<num_skin_particles;p++) {closest_muscle_arr(p) = 0.0;}
    closest_triangle_arr.Resize(num_skin_particles); for(int p=0;p<num_skin_particles;p++) {closest_triangle_arr(p) = 0.0;}
    barycentric_coordinate_arr.Resize(num_skin_particles); for(int p=0;p<num_skin_particles;p++) {barycentric_coordinate_arr(p) = VECTOR_3D<T>(0,0,0);}
    ls_inside.Resize(num_skin_particles); ls_normal.Resize(num_skin_particles);

    //Muscle Levleset
    levelset = new LEVELSET_IMPLICIT_SURFACE<T>(*(new GRID<TV>()), *(new ARRAY<T,VECTOR<int,3> >()));
    char levelset_file[256]; sprintf(levelset_file, "%s/arm_skin_%d.phi", levelset_dir, 0);
    std::istream* levelset_input=FILE_UTILITIES::Safe_Open_Input(levelset_file);levelset->template Read<RW>(*levelset_input);delete levelset_input;
    levelset->Compute_Cell_Minimum_And_Maximum();

    //Muscle Triangles
    Read_Muscle_Tris_From_Quasistatic_Simm_Data(0);

}
//#################################################################
// function Read_Muscle_Tris_From_Quasistatic_Simm_Data
//#################################################################
void Read_Muscle_Tris_From_Quasistatic_Simm_Data(int frame){
    //Muscles
    muscle_tris.Resize(0); 
    DEFORMABLE_OBJECT_LIST_3D<T> *deformable_object_list = new DEFORMABLE_OBJECT_LIST_3D<T>();
    deformable_object_list->template Read_Static_Variables<float>(muscle_directory);
    deformable_object_list->template Read_Dynamic_Variables<float>(muscle_directory, frame);
    int num_musc = deformable_object_list->deformable_objects.m;
    std::cout<<num_musc<<std::endl;
    for(int m=1; m <= num_musc; m++) {
        deformable_object_list->deformable_objects(m)->tetrahedralized_volume->Initialize_Triangulated_Surface();
        muscle_tris.Append(deformable_object_list->deformable_objects(m)->tetrahedralized_volume->triangulated_surface);}
    //Bones
    RIGID_BODY_LIST_3D<T> *rigid_bodies = new RIGID_BODY_LIST_3D<T>(); 
    rigid_bodies->template Read<T>(muscle_directory, frame, true, false);
    RIGID_BODY_LIST_3D<T> *rigid_bodies_1= new RIGID_BODY_LIST_3D<T>(); 
    rigid_bodies_1->template Read<T>(muscle_directory, 0, true, false);
    int start; if (skip_big_ones) start=2; else start=1;
    for(int i=start; i<=rigid_bodies->Number_Of_Elements(); i++) {
        muscle_tris.Append((*rigid_bodies)(i)->triangulated_surface);}
    //Add_Passive_Bones();    //comment this out to skip the fingers and hand
    bool passive;
    for(int i=num_musc+1;i<=muscle_tris.m;i++) {
        int bone_num = i-num_musc;
        if (skip_big_ones) bone_num++;
        if(bone_num > 6) {bone_num=6; passive=true;} else {passive=false;}
        FRAME<T> current_frame((*rigid_bodies)(bone_num)->position, (*rigid_bodies)(bone_num)->orientation);
        FRAME<T>first_frame((*rigid_bodies_1)(bone_num)->position, (*rigid_bodies_1)(bone_num)->orientation);
        for(int p=0;p<muscle_tris(i)->particles.array_collection->Size();p++) {
            if(passive) muscle_tris(i)->particles.X(p)=current_frame*first_frame.Inverse()*muscle_tris(i)->particles.X(p);
            else muscle_tris(i)->particles.X(p)=current_frame*muscle_tris(i)->particles.X(p);}}
    for(int i=0;i<muscle_tris.m;i++) {
        muscle_tris(i)->Refresh_Auxiliary_Structures();
        muscle_tris(i)->Initialize_Triangle_Hierarchy();
        muscle_tris(i)->Update_Bounding_Box();
        muscle_tris(i)->Update_Triangle_List();        }
}

//#####################################################################
// Function Project_Skin_Particles
//#####################################################################
void Project_Skin_Particles(){
    VECTOR_3D<T> old_loc, new_loc;
    for(int i=0;i<num_skin_particles;i++) {
        old_loc = (*((muscle_tris_old(closest_muscle_arr(i)))->triangle_list))(closest_triangle_arr(i)).Point_From_Barycentric_Coordinates(VECTOR_3D<T>(barycentric_coordinate_arr(i)));
        new_loc = (*((muscle_tris(closest_muscle_arr(i)))->triangle_list))(closest_triangle_arr(i)).Point_From_Barycentric_Coordinates(VECTOR_3D<T>(barycentric_coordinate_arr(i)));
        my_skin->particles.X(i) += (new_loc - old_loc);
    }
}

//#####################################################################
// Function Reset_Restlengths, Reset_Velocities and RESET_BENDING 
//#####################################################################
void Reset_Restlengths(){
    (*(solids_parameters.deformable_body_parameters.list(1).linear_springs(1))).Set_Restlength_From_Particles();
    (*(solids_parameters.deformable_body_parameters.list(1).linear_altitude_springs_s3d(1))).Set_Restlength_From_Particles();
}
void Reset_Velocities(){
    TRIANGULATED_SURFACE<T> *my_skin = solids_parameters.deformable_body_parameters.list(1).triangulated_surface;
    for(int c=0;c<my_skin->particles.array_collection->Size();c++) {my_skin->particles.V(c) = VECTOR_3D<T>(0,0,0);}
}
void Reset_Bending(const bool set_rest_angle_to_zero){
    if(set_rest_angle_to_zero) solids_parameters.deformable_body_parameters.list(1).bending_elements(1)->Set_Sine_Half_Rest_Angle(0.0);
    else solids_parameters.deformable_body_parameters.list(1).bending_elements(1)->Set_Constants_From_Particles(skin_bending_stiff,skin_bending_damping);
}
//#####################################################################
// Function Write_Output_Files_Helper
//#####################################################################
void Write_Output_Files_Helper(int f) {
    char levelset_out_file[256];
    Write_Output_Files(f);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/last_frame",f);
    if(write_levelsets){sprintf(levelset_out_file, "Inner_Layer_Sim/output/levelset_%d.phi",f);
    std::ostream *output=FILE_UTILITIES::Safe_Open_Output(levelset_out_file);levelset->template Write<RW>(*output);delete output;}
}

//#####################################################################
};
}
#endif
