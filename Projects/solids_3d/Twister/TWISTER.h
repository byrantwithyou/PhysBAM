//#####################################################################
// Copyright 2002, Robert Bridson, Ronald Fedkiw, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TWISTER
//##################################################################### 
//
// Hanging cloth twisted and untwisted
//
//#####################################################################
// Fedkiw - August 29, 2002
// Bridson - May 21, 2003
// Molino - August 26, 2002
//#####################################################################
#ifndef __TWISTER__
#define __TWISTER__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <fstream>
#include "../CLOTH_EXAMPLE.h"
#include "CUSTOM_EXTERNAL_FORCES_AND_VELOCITIES.h"
#include <Forces_And_Torques/UNIFORM_LINEAR_SPRINGS.h>
namespace PhysBAM{

class TWISTER:public CLOTH_EXAMPLE
{
public:
    int number_side_panels;
    double aspect_ratio,side_length;
    CUSTOM_EXTERNAL_FORCES_AND_VELOCITIES custom_external_forces_and_velocities;

    TWISTER()
                                    :number_side_panels(40),aspect_ratio(2),side_length(1),
                                     custom_external_forces_and_velocities(number_side_panels,aspect_ratio)
    {
        backdoor_id=1;
        final_time=10;
        restart_step_number=0;
        output_directory="Twister/output";
        check_initial_mesh_for_self_intersection=false;
    }

    ~TWISTER()
    {}

//#####################################################################
// Function Initialize_Cloth_State
//#####################################################################
void Initialize_Cloth_State(TRIANGLE_MESH*& triangle_mesh,PARTICLE_3D*& particles,TRIANGULATED_SURFACE*& triangulated_surface,
                                          DEFORMABLE_TRIANGULATED_SURFACE*& cloth)
{
    // shape of cloth mesh
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    double dx=aspect_ratio*side_length/(m-1),dy=side_length/(n-1);
    
    // initialize triangulated surface
    triangle_mesh=new TRIANGLE_MESH();
    particles=new PARTICLE_3D();
    triangulated_surface=new TRIANGULATED_SURFACE(*triangle_mesh,*particles);    
    cloth=new DEFORMABLE_TRIANGULATED_SURFACE(*triangulated_surface);
    if(restart_step_number) Read_Deformable_Triangulated_Surface(*cloth,restart_step_number);
    else{
        triangle_mesh->Initialize_Herring_Bone_Mesh(m,n);
        particles->Update_Position_And_Velocity();particles->Store_Mass();
        for(int k=0;k<triangle_mesh->number_nodes;k++) particles->array_collection->Add_Element();
        double mass_node=aspect_ratio*sqr(side_length)/(m*n);
        copy(mass_node,particles->mass);
        for(int i=0;i<m;i++) for(int j=0;j<n;j++){
            int node=i+m*(j-1);
            particles->X(node)=VECTOR_3D(.5*(i-1)*dx,4*sqr((i-1)/(double)(m-1)-.5),(j-1)*dy);
            particles->V(node)=VECTOR_3D(0,0,0);}}
}
//#####################################################################
// Function Initialize_Cloth_Dynamics
//#####################################################################
void Initialize_Cloth_Dynamics(DEFORMABLE_TRIANGULATED_SURFACE*& cloth)
{
    int k,i,j;
    
    cloth->body_forces.Set_Gravity();
    cloth->Set_External_Forces_And_Velocities(custom_external_forces_and_velocities);
    cloth->Set_CFL_Number(.9);
    cloth->Output_Artificial_Damping_Results();
    
    // shape of cloth mesh
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    
    // make springs
    double spring_length_x=aspect_ratio*side_length/(m-1),spring_length_y=side_length/(n-1);
    
    // stretch springs - x direction
    UNIFORM_LINEAR_SPRINGS *stretch_x=new UNIFORM_LINEAR_SPRINGS(cloth->triangulated_surface.particles);
    stretch_x->segment_mesh.number_nodes=m*n;
    stretch_x->segment_mesh.segments.Resize(2,1,(m-1)*n);
    k=0;
    for(i=0;i<m-1;i++) for(j=0;j<n;j++){
        k++;
        stretch_x->segment_mesh.segments(1,k)=i+m*(j-1);
        stretch_x->segment_mesh.segments(2,k)=i+1+m*(j-1);}
    stretch_x->Set_Stiffness(2);
    stretch_x->Set_Restlength(spring_length_x);
    stretch_x->Set_Overdamping_Fraction(2);
    stretch_x->Artificially_Damp_Max_Strain_Per_Time_Step();
    //stretch_x->Artificially_Damp_Min_And_Max_Strain(true,0,.1);
    cloth->mesh_based_forces.Append(stretch_x);

    // stretch springs - y direction
    UNIFORM_LINEAR_SPRINGS *stretch_y=new UNIFORM_LINEAR_SPRINGS(cloth->triangulated_surface.particles);
    stretch_y->segment_mesh.number_nodes=m*n;
    stretch_y->segment_mesh.segments.Resize(2,1,m*(n-1));
    k=0;
    for(i=0;i<m;i++) for(j=0;j<n-1;j++){
        k++;
        stretch_y->segment_mesh.segments(1,k)=i+m*(j-1);
        stretch_y->segment_mesh.segments(2,k)=i+m*j;}
    stretch_y->Set_Stiffness(2);
    stretch_y->Set_Restlength(spring_length_y);
    stretch_y->Set_Overdamping_Fraction(2); 
    stretch_y->Artificially_Damp_Max_Strain_Per_Time_Step();
    //stretch_y->Artificially_Damp_Min_And_Max_Strain(true,0,.1);
    cloth->mesh_based_forces.Append(stretch_y);
    
    // shear springs
    UNIFORM_LINEAR_SPRINGS *shear=new UNIFORM_LINEAR_SPRINGS(cloth->triangulated_surface.particles);
    shear->segment_mesh.number_nodes=m*n;
    shear->segment_mesh.segments.Resize(2,1,2*(m-1)*(n-1));
    k=0;
    for(i=0;i<m-1;i++) for(j=0;j<n-1;j++){
        k++;
        shear->segment_mesh.segments(1,k)=i+m*(j-1);
        shear->segment_mesh.segments(2,k)=i+1+m*j;
        k++;
        shear->segment_mesh.segments(1,k)=i+1+m*(j-1);
        shear->segment_mesh.segments(2,k)=i+m*j;}
    shear->Set_Stiffness(2/sqrt(2.));
    shear->Set_Restlength(sqrt(sqr(spring_length_x)+sqr(spring_length_y)));
    shear->Set_Overdamping_Fraction(2); 
    shear->Artificially_Damp_Max_Strain_Per_Time_Step();
    //shear->Artificially_Damp_Min_And_Max_Strain(true,0,.1);
    cloth->mesh_based_forces.Append(shear);
    
    // triangle_bending_springs
    TRIANGLE_BENDING_SPRINGS* bend=new TRIANGLE_BENDING_SPRINGS(cloth->triangulated_surface.particles);
    bend->Set_Quadruples_From_Triangle_Mesh(cloth->triangulated_surface.triangle_mesh);
    bend->Set_Stiffness(.002);
    bend->Set_Damping(.002);
    cloth->mesh_based_forces.Append(bend);

    // set up repulsion springs for collisions !!!!!!
    repulsion_springs_initialized=true;
    collisions_repulsion_spring_youngs_modulus=stretch_x->youngs_modulus;
    collisions_repulsion_spring_restlength=min(spring_length_x,spring_length_y);
    double mass_node=aspect_ratio*sqr(side_length)/(m*n);
    collisions_repulsion_spring_mass=mass_node;
}
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies()
{
    int index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/medium_cylinder",(T).2);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->position=VECTOR_3D(.5,.15,.5);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->orientation=QUATERNION(pi/2,VECTOR_3D(1,0,0));
    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Update_Object_Positions
//#####################################################################
void Update_Object_Positions(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const double time)
{
    rigid_bodies(1)->position=VECTOR_3D(.5,.15,.5);
    rigid_bodies(1)->orientation=QUATERNION(time,VECTOR_3D(0,1,0))*QUATERNION(pi/2,VECTOR_3D(1,0,0));
}
//#####################################################################
// Function Update_Object_Velocities
//#####################################################################
void Update_Object_Velocities(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const double time)
{
    rigid_bodies(1)->angular_velocity=VECTOR_3D(0,1,0);
}
//#####################################################################
};
}
#endif

