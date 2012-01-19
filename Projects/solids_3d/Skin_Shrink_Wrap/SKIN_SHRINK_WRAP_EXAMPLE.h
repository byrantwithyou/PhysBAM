//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Unnur Gretarsdottir.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHRINK_WRAP_EXAMPLE
//#####################################################################
#ifndef __SKIN_SHRINK_WRAP_EXAMPLE__
#define __SKIN_SHRINK_WRAP_EXAMPLE__

//Fix las shrink
//Fix Shrink to bone axis
//constrain to normal of levelset is ok?? 
//fix force to levelset
//postprocess substep is in wrong spot & output at each frame is hacked  

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>

namespace PhysBAM{
template<class T,class RW>
class SKIN_SHRINK_WRAP_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    std::string input_file;std::string inner_layer_file;
    int substep_counter;
    int shrink_every_x_steps;
    T shrink_percentage;
    T force_to_levelset;
    int grow_steps, shrink_steps, settle_steps;
    int mode; //0-growing  1-shrinking  2-settling
    ARRAY<bool> ever_inside;
    ARRAY<VECTOR_3D<T> > ls_normal;
    T skin_ls_stiff, skin_ls_damping, skin_las_stiff, skin_las_damping, skin_bending_stiff, skin_bending_damping, bending_cutoff_ratio;

    SKIN_SHRINK_WRAP_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE), substep_counter(0), mode(0)
    {   
        /////////////////////Tuning Parameters/////////////////////////
        shrink_every_x_steps = 20; //In Shrinking mode, we call Shrink(), every X steps
        shrink_percentage = -.2;   //When restlenghts shrink, they reduce their length by this number (should be negative)
        force_to_levelset = .02;   //When settling, we apply a slight force to the levelset for extra definition
        //number of substeps to allow mesh to grow, shrink and settle      
        grow_steps = 200;    
        shrink_steps = 0;        
        settle_steps = 100;        
        //stiffness parameters - note the skin density is currently 1 - probably not good...
        skin_ls_stiff = .2/(1+sqrt((T)2)); skin_ls_damping = 2;
        skin_las_stiff = .2*4/(1+sqrt((T)2)); skin_las_damping = 4;
        skin_bending_stiff = 1e-3; skin_bending_damping = 1e-3;
        bending_cutoff_ratio = .5;
        /////////////////////Other Parameters/////////////////////////
        frame_rate =.1; last_frame=1; restart=false; //we output manually in Postprocess_Substep, so frames are pretty irrelevant
        write_output_files = false;
        solids_parameters.use_constant_mass = true;
        solids_parameters.cfl=(T)5.9;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        solids_parameters.perform_self_collision=false;
        solids_parameters.collide_with_interior=false;
        output_directory="Skin_Shrink_Wrap/output";
        inner_layer_file="Skin_Shrink_Wrap/inner_layer.tri";
        input_file="Skin_Shrink_Wrap/arm_skin.tri.gz";
      }
    ~SKIN_SHRINK_WRAP_EXAMPLE()
    {}

//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
//this is called incorrectly in the Physbam version of Solids_Evolution_3D.cpp - make sure to fix it!!!!
void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    Update_Ever_Inside(); substep_counter++;
    if(mode == 0) {Reset_Restlengths();std::cout<<"GROWING: ";} //reset restlengths as we grow
    if(mode == 1) {if(((substep_counter-1) % shrink_every_x_steps) == 0) Shrink();std::cout<<"SHRINKING: ";} //Shrink if it's the right time
    if(mode == 2) {std::cout<<"SETTLING: ";} //settling - do nothing except apply a slight force to the levelset (done in Add_External_Forces
    if(substep_counter == 1){mode = 0; Reset_Restlengths();Reset_Velocities();Reset_Ever_Inside();Set_Bending();}//switch to mode 1
    if(substep_counter == grow_steps){mode = 1; Reset_Restlengths();Reset_Velocities();Reset_Ever_Inside();}    //switch to mode 2
    if(substep_counter == shrink_steps+grow_steps) {mode =2; Reset_Restlengths();Reset_Velocities();Reset_Ever_Inside();} //switch to mode 3
    if(substep_counter == shrink_steps+grow_steps+settle_steps) {Finished();} //finished - write out tri file and quit
    std::cout<<"Substep number "<<substep_counter-1<<" completed"<<std::endl;
    Write_Output_Files(substep_counter); FILE_UTILITIES::Write_To_Text_File(output_directory+"/last_frame",substep_counter);
}
//#####################################################################
// Function Finished
//#####################################################################
void Finished(){
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(inner_layer_file);
    solids_parameters.deformable_body_parameters.list(1).triangulated_surface->template Write<RW>(*output);delete output;
    exit(0);
}
//#####################################################################
// Function Reset_Restlengths and Reset_Velocities 
//#####################################################################
void Reset_Restlengths(){
    (*(solids_parameters.deformable_body_parameters.list(1).linear_springs(1))).Set_Restlength_From_Particles();
    (*(solids_parameters.deformable_body_parameters.list(1).linear_altitude_springs_s3d(1))).Set_Restlength_From_Particles();}
void Reset_Velocities(){
    TRIANGULATED_SURFACE<T> *my_skin = solids_parameters.deformable_body_parameters.list(1).triangulated_surface;
    for(int c=0;c<my_skin->particles.array_collection->Size();c++) {my_skin->particles.V(c) = VECTOR_3D<T>(0,0,0);}}
void Set_Bending(){
    (*(solids_parameters.deformable_body_parameters.list(1).bending_elements(1))).Set_Sine_Half_Rest_Angle(0.0);}

//#####################################################################
// Function Ever_Inside Stuff 
//#####################################################################
void Reset_Ever_Inside(){
    for(int i=0;i<ever_inside.m;i++) {ever_inside(i) = false; ls_normal(i) = VECTOR_3D<T>(0.0,0.0,0.0);}}
void Update_Ever_Inside() {
    std::cout<<"Updating Ever Inside"<<std::endl;
    ARRAY<bool> is_inside = solids_parameters.deformable_body_parameters.list(1).collisions.enforce_collision;
    for(int i=0;i<is_inside.m;i++) {
        if(is_inside(i)) {
            ever_inside(i) = true;
            ls_normal(i) = solids_parameters.deformable_body_parameters.list(1).collisions.collision_normal(i);}}
}
//#####################################################################
// Shrink
// ===============
// Shrinks restlengths - note that if both nodes on a segment are now
// touching the levelset, we don't want to shrink it anymore, just 
// leave it be (so set it's restlenght to the current length). You could
// change this and keep shrinking these segments, or not reset their
// restlengths, but this approach worked best for me.
//#####################################################################
void Shrink(){
    DEFORMABLE_OBJECT_3D<T> *my_def_obj = &(solids_parameters.deformable_body_parameters.list(1));
    LINEAR_SPRINGS<T,VECTOR_3D<T> > *ls = my_def_obj->linear_springs(1);
    LINEAR_ALTITUDE_SPRINGS_3D<T> *las = my_def_obj->linear_altitude_springs_3d(1);
    for(int u=0;u<ls->restlength.m;u++) {
        int node1=ls->segment_mesh.segments(1,u);
        int node2=ls->segment_mesh.segments(2,u);
        if(ever_inside(node1)&&ever_inside(node2)) { //both nodes are on the levelset
            ls->restlength(u) = (my_def_obj->particles.X(node1) - my_def_obj->particles.X(node2)).Magnitude();continue;}
        VECTOR_3D<T> spring = (my_def_obj->particles.X(node1) - my_def_obj->particles.X(node2)).Normalized();
        //example of anisotropic shrinking (this is along the z-axis, but you probably don't want this for these sims
        //VECTOR_3D<T> axis(0.0,0.0,1.0);
        //T ratio = (spring.Projected_On_Unit_Direction(axis)).Magnitude();
        //ratio *= ratio; ratio = 1.0-ratio;
        T ratio = 1.0;  //Don't want to shrink anisotropically in these sims, so just set ratio to 1
        ratio = 1.0 + (ratio * shrink_percentage); ls->restlength(u) *= ratio;  }
}
//#####################################################################
// Function Set_External_Velocities and Zero_Out_Enslaved_Velocity_Nodes  - why does this compile?????? no ls_inside defined
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time) {
    if(mode > 1) {   //Constrain direction normal to levelset to be 0
        for(int c=0;c<V.m;c++) {
            if(ever_inside(c)) { //nodes on levelset stay there...
                V(c)-=VECTOR_3D<T>::Dot_Product(V(c),ls_normal(c))*ls_normal(c); 
            }}
}}
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time) {    
    if(mode > 1) {   //Constrain direction normal to levelset to be 0
        for(int c=0;c<V.m;c++) {
            if(ever_inside(c)) { 
                V(c)-=VECTOR_3D<T>::Dot_Product(V(c),ls_normal(c))*ls_normal(c); 
            }}
}}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
void Add_External_Forces(ARRAY<VECTOR_3D<T> >& F,const T time) PHYSBAM_OVERRIDE {
    if(mode!= 2) return;
    IMPLICIT_SURFACE<T> *levelset1 = solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->implicit_surface;
    for(int p=0;p<F.m;p++) {
        VECTOR_3D<T> current_location = solids_parameters.deformable_body_parameters.list(1).triangulated_surface->particles.X(p);
        VECTOR_3D<T> n = levelset1->Extended_Normal(current_location);
        T depth; levelset1->Lazy_Outside_Extended_Levelset_And_Value(current_location, depth); 
        levelset1->Lazy_Inside_Extended_Levelset_And_Value(current_location, depth); 
        T factor = force_to_levelset;  //this could depend on distance from levelset, bending strength or anything else - for now it's just a constant force
        F(p) += (-n*factor*depth);
    }
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data() {
    //Skin Mesh
     SKIN_DEFORMABLE_OBJECT_3D<T> *my_skin = new SKIN_DEFORMABLE_OBJECT_3D<T>(); my_skin->example = this;
    solids_parameters.deformable_body_parameters.list.deformable_objects.Append(my_skin);
    int index = solids_parameters.deformable_body_parameters.list.deformable_objects.m;
    solids_parameters.deformable_body_parameters.list.deformable_objects(index)->Allocate_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).triangulated_surface;
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_file);triangulated_surface.template Read<RW>(*input);delete input;
    std::cout << "total vertices = " << particles.array_size << std::endl;std::cout << "total triangles = " << triangle_mesh.triangles.m << std::endl;
    particles.store_velocity = false; particles.Update_Velocity();  // forcing velocity resize - not sure why this is needed
    particles.store_mass = false; particles.Store_Mass(); // forcing mass resize - not sure why
    for(int p=0;p<particles.array_collection->Size();p++) {particles.mass(p) = 0.0;}
    triangulated_surface.Set_Density(1);triangulated_surface.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    triangulated_surface.Update_Bounding_Box();
    ever_inside.Resize(solids_parameters.deformable_body_parameters.list(index).triangulated_surface->particles.array_collection->Size());
    ls_normal.Resize(solids_parameters.deformable_body_parameters.list(index).triangulated_surface->particles.array_collection->Size());
    //Muscle Levleset
    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("Skin_Shrink_Wrap/Arm_Output/arm_skin", 1, true, true, true);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() {
    Get_Initial_Data();
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(1).triangulated_surface;
//    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface);
    solids_parameters.deformable_body_parameters.list(1).Add_Edge_Springs(triangulated_surface.triangle_mesh,skin_ls_stiff,skin_ls_damping);
    solids_parameters.deformable_body_parameters.list(1).Add_Altitude_Springs(triangulated_surface.triangle_mesh,skin_las_stiff,skin_las_damping);
    solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements(triangulated_surface.triangle_mesh,skin_bending_stiff,skin_bending_damping);
    solids_parameters.deformable_body_parameters.list(1).bending_elements(1)->Set_Area_Cutoff_From_Triangulated_Surface(triangulated_surface,bending_cutoff_ratio);
}
//#####################################################################
};
}
#endif
