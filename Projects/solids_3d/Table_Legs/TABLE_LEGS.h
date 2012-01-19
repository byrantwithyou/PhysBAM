//#####################################################################
// Copyright 2002, Robert Bridson, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TABLE_LEGS
//##################################################################### 
//
// Cloth draped over 4 table legs, a ball drags it down to the floor
//
//#####################################################################
// Fedkiw - June 26, 2002
// Bridson - June 28, 2002
//#####################################################################
#ifndef __TABLE_LEGS__
#define __TABLE_LEGS__

#include <fstream>
#include "CLOTH_EXAMPLE.h"
#include <Forces_And_Torques/UNIFORM_BENDING_SPRINGS.h>
#include <Forces_And_Torques/UNIFORM_LINEAR_SPRINGS.h>
#include <Rigid_Bodies/RIGID_PLANE.h>
#include <Rigid_Bodies/RIGID_SPHERE.h>
#include <Rigid_Bodies/RIGID_TRIANGULATED_SOLID.h>

namespace PhysBAM {
class TABLE_LEGS: public CLOTH_EXAMPLE
{
public:
    int number_side_panels,m,n;
private:
    TRIANGLE_MESH leg_mesh;
    PARTICLE_3D leg1_particles,leg2_particles,leg3_particles,leg4_particles;
    TRIANGULATED_SURFACE leg1_surf,leg2_surf,leg3_surf,leg4_surf;

public:
    TABLE_LEGS(int number_side_panels_input)
        :number_side_panels(number_side_panels_input),leg1_surf(leg_mesh,leg1_particles),leg2_surf(leg_mesh,leg2_particles),leg3_surf(leg_mesh,leg3_particles),leg4_surf(leg_mesh,leg4_particles)
    {backdoor_id=2;}

    ~TABLE_LEGS()
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
    double aspect_ratio=1,side_length=2;
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
        for(int i=0;i<m;i++) for(int j=0;j<n;j++){
            int node=i+m*(j-1);
            particles->X(node)=VECTOR_3D((i-1)*dx-.5*aspect_ratio*side_length,1.1,(j-1)*dy-.5*side_length);
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
    triangle_collisions->Set_Up_Repulsion_Spring(stretch_x->youngs_modulus,min(spring_length_x,spring_length_y),mass_node);
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
    scripted_objects.Resize(1,6);
    scripted_objects_friction.Resize(1,6);
    // plane
    scripted_objects(1)=new RIGID_PLANE();
    scripted_objects_friction(1)=.1;
    // table leg mesh
    leg_mesh.triangles.Exact_Resize(3,1,20);
    leg_mesh.triangles(1,1)=1;  leg_mesh.triangles(2,1)=2;  leg_mesh.triangles(3,1)=3;
    leg_mesh.triangles(1,2)=1;  leg_mesh.triangles(2,2)=3;  leg_mesh.triangles(3,2)=4;
    leg_mesh.triangles(1,3)=1;  leg_mesh.triangles(2,3)=5;  leg_mesh.triangles(3,3)=2;
    leg_mesh.triangles(1,4)=5;  leg_mesh.triangles(2,4)=6;  leg_mesh.triangles(3,4)=2;
    leg_mesh.triangles(1,5)=2;  leg_mesh.triangles(2,5)=6;  leg_mesh.triangles(3,5)=3;
    leg_mesh.triangles(1,6)=6;  leg_mesh.triangles(2,6)=7;  leg_mesh.triangles(3,6)=3;
    leg_mesh.triangles(1,7)=3;  leg_mesh.triangles(2,7)=7;  leg_mesh.triangles(3,7)=4;
    leg_mesh.triangles(1,8)=7;  leg_mesh.triangles(2,8)=8;  leg_mesh.triangles(3,8)=4;
    leg_mesh.triangles(1,9)=4;  leg_mesh.triangles(2,9)=8;  leg_mesh.triangles(3,9)=1;
    leg_mesh.triangles(1,10)=8; leg_mesh.triangles(2,10)=5; leg_mesh.triangles(3,10)=1;
    leg_mesh.triangles(1,11)=5; leg_mesh.triangles(2,11)=9; leg_mesh.triangles(3,11)=6;
    leg_mesh.triangles(1,12)=9; leg_mesh.triangles(2,12)=10;leg_mesh.triangles(3,12)=6;
    leg_mesh.triangles(1,13)=6; leg_mesh.triangles(2,13)=10;leg_mesh.triangles(3,13)=7;
    leg_mesh.triangles(1,14)=10;leg_mesh.triangles(2,14)=11;leg_mesh.triangles(3,14)=7;
    leg_mesh.triangles(1,15)=7; leg_mesh.triangles(2,15)=11;leg_mesh.triangles(3,15)=8;
    leg_mesh.triangles(1,16)=11;leg_mesh.triangles(2,16)=12;leg_mesh.triangles(3,16)=8;
    leg_mesh.triangles(1,17)=8; leg_mesh.triangles(2,17)=12;leg_mesh.triangles(3,17)=5;
    leg_mesh.triangles(1,18)=12;leg_mesh.triangles(2,18)=9; leg_mesh.triangles(3,18)=5;
    leg_mesh.triangles(1,19)=9; leg_mesh.triangles(2,19)=11;leg_mesh.triangles(3,19)=10;
    leg_mesh.triangles(1,20)=9; leg_mesh.triangles(2,20)=12;leg_mesh.triangles(3,20)=11;
    leg_mesh.number_nodes=12;
    // first table leg
    leg1_particles.Store_Position();
    leg1_particles.Set_Array_Buffer_Size(12);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(1)=VECTOR_3D(.25,-2,.25);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(2)=VECTOR_3D(.4,-2,.25);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(3)=VECTOR_3D(.4,-2,.4);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(4)=VECTOR_3D(.25,-2,.4);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(5)=VECTOR_3D(.25,.965,.25);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(6)=VECTOR_3D(.4,.965,.25);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(7)=VECTOR_3D(.4,.965,.4);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(8)=VECTOR_3D(.25,.965,.4);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(9)=VECTOR_3D(.28,1,.28);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(10)=VECTOR_3D(.37,1,.28);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(11)=VECTOR_3D(.37,1,.37);
    leg1_particles.array_collection->Add_Element(); leg1_particles.X(12)=VECTOR_3D(.28,1,.37);
    scripted_objects(2)=new RIGID_TRIANGULATED_SOLID(leg1_surf);
    scripted_objects_friction(2)=0;
    // second table leg
    leg2_particles.Store_Position();
    leg2_particles.Set_Array_Buffer_Size(12);
    for(int i=0;i<12;i++) {leg2_particles.array_collection->Add_Element(); leg2_particles.X(i)=leg1_particles.X(i)+VECTOR_3D(-.65,0,0);}
    scripted_objects(3)=new RIGID_TRIANGULATED_SOLID(leg2_surf);
    scripted_objects_friction(3)=0;
    // third table leg
    leg3_particles.Store_Position();
    leg3_particles.Set_Array_Buffer_Size(12);
    for(int i=0;i<12;i++) {leg3_particles.array_collection->Add_Element(); leg3_particles.X(i)=leg2_particles.X(i)+VECTOR_3D(0,0,-.65);}
    scripted_objects(4)=new RIGID_TRIANGULATED_SOLID(leg3_surf);
    scripted_objects_friction(4)=0;
    // fourth table leg
    leg4_particles.Store_Position();
    leg4_particles.Set_Array_Buffer_Size(12);
    for(int i=0;i<12;i++) {leg4_particles.array_collection->Add_Element(); leg4_particles.X(i)=leg3_particles.X(i)+VECTOR_3D(.65,0,0);}
    scripted_objects(5)=new RIGID_TRIANGULATED_SOLID(leg4_surf);
    scripted_objects_friction(5)=0;
    // sphere
    scripted_objects(6)=new RIGID_SPHERE(VECTOR_3D(0,3,0),.2);
    scripted_objects_friction(6)=.3;
    // since we're scripting the objects, make inertia tensor identity to make angular_velocity=angular_momentum
    for(int i=1;i<=scripted_objects.m;++i)scripted_objects(i)->moments_of_inertia=VECTOR_3D(1,1,1);
}
//#####################################################################
// Function Update_Object_Positions
//#####################################################################
void Update_Object_Positions(ARRAY<RIGID_BODY<TV>*>& scripted_objects,double time)
{
    if(time <= 2.79) scripted_objects(6)->position=VECTOR_3D(0,3-time,0);
    else scripted_objects(6)->position=VECTOR_3D(time-2.79,.21,0);
}
//#####################################################################
// Function Update_Object_Velocities
//#####################################################################
void Update_Object_Velocities(ARRAY<RIGID_BODY<TV>*>& scripted_objects,double time)
{
    if(time <= 2.79) scripted_objects(6)->velocity=VECTOR_3D(0,-1,0);
    else scripted_objects(6)->velocity=VECTOR_3D(1,0,0);
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
