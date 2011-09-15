//#####################################################################
// Copyright 2002, Robert Bridson, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPINNING_BALL
//##################################################################### 
//
// Cloth draped over a spinning ball and a bump
//
//#####################################################################
// Fedkiw - June 26, 2002
// Bridson - June 28, 2002
//#####################################################################
#ifndef __SPINNING_BALL__
#define __SPINNING_BALL__

#include <fstream>
#include "CLOTH_EXAMPLE.h"
#include <Forces_And_Torques/UNIFORM_BENDING_SPRINGS.h>
#include <Forces_And_Torques/UNIFORM_LINEAR_SPRINGS.h>
#include <Rigid_Bodies/RIGID_PLANE.h>
#include <Rigid_Bodies/RIGID_SPHERE.h>
#include <Rigid_Bodies/RIGID_TRIANGULATED_SOLID.h>

namespace PhysBAM {
class SPINNING_BALL: public CLOTH_EXAMPLE
{
public:
    int number_side_panels,m,n;
private:
    TRIANGLE_MESH bump_mesh;
    PARTICLE_3D bump_particles;
    TRIANGULATED_SURFACE bump_surf;

public:
    SPINNING_BALL(int number_side_panels_input)
        :number_side_panels(number_side_panels_input),bump_surf(bump_mesh,bump_particles)
    {backdoor_id=2;}

    ~SPINNING_BALL()
    {}

//#####################################################################
// Function Initialize_Time_Stepping
//#####################################################################
void Initialize_Time_Stepping(double &cfl,double &final_time)
{
    cfl=.9;final_time=9;
}
//#####################################################################
// Function Initialize_Cloth
//#####################################################################
void Initialize_Cloth(TRIANGLE_MESH*& triangle_mesh,PARTICLE_3D*& particles,TRIANGULATED_SURFACE*& triangulated_surface,DEFORMABLE_TRIANGULATED_SURFACE*& cloth,TRIANGLE_COLLISIONS*& triangle_collisions,int restart_step_number,const char* data_directory)
{
    int i,j,k;

    // shape of cloth mesh
    double aspect_ratio=1,side_length=1.6;
    m=(int)(aspect_ratio*number_side_panels)+1;n=number_side_panels+1;
    double dx=aspect_ratio*side_length/(m-1),dy=side_length/(n-1);
    
    // initialize triangulated surface
    triangle_mesh=new TRIANGLE_MESH();
    particles=new PARTICLE_3D();
    triangulated_surface=new TRIANGULATED_SURFACE(*triangle_mesh,*particles);
    triangle_mesh->Initialize_Square_Mesh(m,n);
    particles->Update_Position_And_Velocity();particles->Store_Mass();
    for(k=1;k<=triangle_mesh->number_nodes;k++) particles->array_collection->Add_Element();
    if(!restart_step_number){    
        for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
            int node=i+m*(j-1);
            particles->X(node)=VECTOR_3D((i-1)*dx-.5*aspect_ratio*side_length,.5,(j-1)*dy-.5*side_length);
            particles->V(node)=VECTOR_3D(0,0,0);}}
    else{
        std::fstream input;
        char File_Name[256];sprintf(File_Name,"%s/position.%d",data_directory,restart_step_number);
        input.open(File_Name,std::ios::in|std::ios::binary);
        double data_double;
        for(k=1;k<=triangle_mesh->number_nodes;k++){
            input.read((char*)&data_double,8);particles->X(k).x=(double)data_double;
            input.read((char*)&data_double,8);particles->X(k).y=(double)data_double;
            input.read((char*)&data_double,8);particles->X(k).z=(double)data_double;}
        input.close();
        sprintf(File_Name,"%s/velocity.%d",data_directory,restart_step_number);
        input.open(File_Name,std::ios::in|std::ios::binary);
        for(k=1;k<=triangle_mesh->number_nodes;k++){
            input.read((char*)&data_double,8);particles->V(k).x=(double)data_double;
            input.read((char*)&data_double,8);particles->V(k).y=(double)data_double;
            input.read((char*)&data_double,8);particles->V(k).z=(double)data_double;}
        input.close();}
    // set mass of the particles - all equal mass
    double mass_node=aspect_ratio*sqr(side_length)/(m*n);
    copy(mass_node,particles->mass);

    cloth=new DEFORMABLE_TRIANGULATED_SURFACE(*triangulated_surface);
    cloth->Set_Thickness(.001);
    // make springs
    double spring_length_x=aspect_ratio*side_length/(m-1),spring_length_y=side_length/(n-1);
    // stretch springs - x direction
    UNIFORM_LINEAR_SPRINGS *stretch_x=new UNIFORM_LINEAR_SPRINGS(*particles);
    stretch_x->segment_mesh.number_nodes=m*n;
    stretch_x->segment_mesh.segments.Resize(2,1,(m-1)*n);
    k=0;
    for(i=1;i<=m-1;i++) for(j=1;j<=n;j++){
        k++;
        stretch_x->segment_mesh.segments(1,k)=i+m*(j-1);
        stretch_x->segment_mesh.segments(2,k)=i+1+m*(j-1);}
    stretch_x->Set_Stiffness(2*1/(1+sqrt(2.)));
    stretch_x->Set_Restlength(spring_length_x);
    stretch_x->Set_Overdamping_Fraction(2); 
    stretch_x->Set_Artificial_Max_Strain_Per_Time_Step(.1);
    stretch_x->Set_Artificial_Min_Strain(0);
    stretch_x->Set_Artificial_Max_Strain(.1);
    cloth->mesh_based_forces.Append(stretch_x);
    // stretch springs - y direction
    UNIFORM_LINEAR_SPRINGS *stretch_y=new UNIFORM_LINEAR_SPRINGS(*particles);
    stretch_y->segment_mesh.number_nodes=m*n;
    stretch_y->segment_mesh.segments.Resize(2,1,m*(n-1));
    k=0;
    for(i=1;i<=m;i++) for(j=1;j<=n-1;j++){
        k++;
        stretch_y->segment_mesh.segments(1,k)=i+m*(j-1);
        stretch_y->segment_mesh.segments(2,k)=i+m*j;}
    stretch_y->Set_Stiffness(2*1/(1+sqrt(2.)));
    stretch_y->Set_Restlength(spring_length_y);
    stretch_y->Set_Overdamping_Fraction(2); 
    stretch_y->Set_Artificial_Max_Strain_Per_Time_Step(.1);
    stretch_y->Set_Artificial_Min_Strain(0);
    stretch_y->Set_Artificial_Max_Strain(.1);
    cloth->mesh_based_forces.Append(stretch_y);
    // shear springs
    UNIFORM_LINEAR_SPRINGS *shear=new UNIFORM_LINEAR_SPRINGS(*particles);
    shear->segment_mesh.number_nodes=m*n;
    shear->segment_mesh.segments.Resize(2,1,2*(m-1)*(n-1));
    k=0;
    for(i=1;i<=m-1;i++) for(j=1;j<=n-1;j++){
        k++;
        shear->segment_mesh.segments(1,k)=i+m*(j-1);
        shear->segment_mesh.segments(2,k)=i+1+m*j;
        k++;
        shear->segment_mesh.segments(1,k)=i+1+m*(j-1);
        shear->segment_mesh.segments(2,k)=i+m*j;}
    shear->Set_Stiffness(2*1/(1+sqrt(2.)));
    shear->Set_Restlength(sqrt(sqr(spring_length_x)+sqr(spring_length_y)));
    shear->Set_Overdamping_Fraction(2); 
    shear->Set_Artificial_Max_Strain_Per_Time_Step(.1);
    shear->Set_Artificial_Min_Strain(0);
    shear->Set_Artificial_Max_Strain(.1);
    cloth->mesh_based_forces.Append(shear);
    // bending springs - x direction
    UNIFORM_BENDING_SPRINGS *bend_x=new UNIFORM_BENDING_SPRINGS(*particles);
    bend_x->bending_triples.Resize(3,1,(m-2)*n);
    k=0;
    for(i=2;i<=m-1;i++) for(j=1;j<=n;j++){
        k++;
        bend_x->bending_triples(1,k)=i-1+m*(j-1);
        bend_x->bending_triples(2,k)=i+m*(j-1);
        bend_x->bending_triples(3,k)=i+1+m*(j-1);}
    bend_x->Set_Stiffness(1./40);
    bend_x->Set_Restlength(2*spring_length_x);
    bend_x->Set_Overdamping_Fraction(4); 
    cloth->mesh_based_forces.Append(bend_x);
    // bending springs - y direction
    UNIFORM_BENDING_SPRINGS *bend_y=new UNIFORM_BENDING_SPRINGS(*particles);
    bend_y->bending_triples.Resize(3,1,m*(n-2));
    k=0;
    for(i=1;i<=m;i++) for(j=2;j<=n-1;j++){
        k++;
        bend_y->bending_triples(1,k)=i+m*(j-2);
        bend_y->bending_triples(2,k)=i+m*(j-1);
        bend_y->bending_triples(3,k)=i+m*j;}
    bend_y->Set_Stiffness(1./40);
    bend_y->Set_Restlength(2*spring_length_y);
    bend_y->Set_Overdamping_Fraction(4); 
    cloth->mesh_based_forces.Append(bend_y);

    // triangle collisions and friction - note that the initial bounding boxes are initialized
    triangle_collisions=new TRIANGLE_COLLISIONS(*triangulated_surface);
    triangle_collisions->Set_Small_Number(1e-8);
    triangle_collisions->Set_Repulsion_Thickness(cloth->thickness);
    triangle_collisions->Set_Collision_Thickness(1e-6);
    triangle_collisions->Set_Up_Repulsion_Spring(stretch_x->youngs_modulus/min(spring_length_x,spring_length_y),mass_node);
    triangle_collisions->Output_Repulsion_Results();
    triangle_collisions->Output_Collision_Results();
    //triangle_collisions->Output_Number_Checked();
    triangle_collisions->Set_Friction_Coefficient(.4);
    triangle_collisions->Set_Attempts_For_Nonrigid_Collisions(10);
}
//#####################################################################
// Function Initialize_Scripted_Objects
//#####################################################################
void Initialize_Scripted_Objects(ARRAY<RIGID_BODY<TV>*>& scripted_objects,ARRAY<double>& scripted_objects_friction)
{
    scripted_objects.Resize(3);
    scripted_objects_friction.Resize(3);
    // plane
    scripted_objects(1)=new RIGID_PLANE();
    scripted_objects_friction(1)=.15;
    // bump
    bump_mesh.triangles.Exact_Resize(3,1,12);
    bump_mesh.triangles(1,1)=1; bump_mesh.triangles(2,1)=2; bump_mesh.triangles(3,1)=3;
    bump_mesh.triangles(1,2)=1; bump_mesh.triangles(2,2)=3; bump_mesh.triangles(3,2)=4;
    bump_mesh.triangles(1,3)=5; bump_mesh.triangles(2,3)=7; bump_mesh.triangles(3,3)=6;
    bump_mesh.triangles(1,4)=5; bump_mesh.triangles(2,4)=8; bump_mesh.triangles(3,4)=7;
    bump_mesh.triangles(1,5)=1; bump_mesh.triangles(2,5)=5; bump_mesh.triangles(3,5)=2;
    bump_mesh.triangles(1,6)=5; bump_mesh.triangles(2,6)=6; bump_mesh.triangles(3,6)=2;
    bump_mesh.triangles(1,7)=2; bump_mesh.triangles(2,7)=6; bump_mesh.triangles(3,7)=3;
    bump_mesh.triangles(1,8)=6; bump_mesh.triangles(2,8)=7; bump_mesh.triangles(3,8)=3;
    bump_mesh.triangles(1,9)=3; bump_mesh.triangles(2,9)=7; bump_mesh.triangles(3,9)=4;
    bump_mesh.triangles(1,10)=7;bump_mesh.triangles(2,10)=8;bump_mesh.triangles(3,10)=4;
    bump_mesh.triangles(1,11)=4;bump_mesh.triangles(2,11)=8;bump_mesh.triangles(3,11)=1;
    bump_mesh.triangles(1,12)=8;bump_mesh.triangles(2,12)=5;bump_mesh.triangles(3,12)=1;
    bump_mesh.number_nodes=8;
    bump_particles.Store_Position();
    bump_particles.Set_Array_Buffer_Size(8);
    bump_particles.array_collection->Add_Element(); bump_particles.X(1)=VECTOR_3D(.1,-.3,.1);
    bump_particles.array_collection->Add_Element(); bump_particles.X(2)=VECTOR_3D(.9,-.3,.1);
    bump_particles.array_collection->Add_Element(); bump_particles.X(3)=VECTOR_3D(.9,-.3,.9);
    bump_particles.array_collection->Add_Element(); bump_particles.X(4)=VECTOR_3D(.1,-.3,.9);
    bump_particles.array_collection->Add_Element(); bump_particles.X(5)=VECTOR_3D(.4,.2,.4);
    bump_particles.array_collection->Add_Element(); bump_particles.X(6)=VECTOR_3D(.6,.2,.4);
    bump_particles.array_collection->Add_Element(); bump_particles.X(7)=VECTOR_3D(.6,.2,.6);
    bump_particles.array_collection->Add_Element(); bump_particles.X(8)=VECTOR_3D(.4,.2,.6);
    scripted_objects(2)=new RIGID_TRIANGULATED_SOLID(bump_surf);
    scripted_objects_friction(2)=.2;
    // rotating sphere
    scripted_objects(3)=new RIGID_SPHERE(VECTOR_3D(0,.25,0),.2);
    scripted_objects_friction(3)=.5;
}
//#####################################################################
// Function Update_Object_Positions
//#####################################################################
void Update_Object_Positions(ARRAY<RIGID_BODY<TV>*>& scripted_objects,double time)
{
}
//#####################################################################
// Function Update_Object_Velocities
//#####################################################################
void Update_Object_Velocities(ARRAY<RIGID_BODY<TV>*>& scripted_objects,double time)
{
    scripted_objects(3)->angular_velocity=VECTOR_3D(0,min(2,.5*time),0);
}
//#####################################################################
// Function Set_External_Forces
//#####################################################################
void Set_External_Forces(ARRAY<VECTOR_3D>& F,ARRAY<double>& mass,double time)
{
    int k;
    for(k=1;k<=F.m;k++) F(k)=VECTOR_3D(0,0,0); // initialize
    for(k=1;k<=F.m;k++) F(k).y+=(-9.8*mass(k)); // add gravity
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D>& V,double time,bool inhomogeneous=true)
{
}
//#####################################################################
};
}
#endif
