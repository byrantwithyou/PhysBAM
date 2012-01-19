//#####################################################################
// Copyright 2003, Robert Bridson.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_BALL
//##################################################################### 
//
// Triangular cloth draped over a ball (for convergence tests)
//
//#####################################################################
// Bridson - February 13, 2003
//#####################################################################
#ifndef __SIMPLE_BALL__
#define __SIMPLE_BALL__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <fstream>
#include "../CLOTH_EXAMPLE.h"
#include <Forces_And_Torques/UNIFORM_LINEAR_SPRINGS.h>
namespace PhysBAM{

class SIMPLE_BALL:public CLOTH_EXAMPLE
{
public:
    int number_subdivisions;

    SIMPLE_BALL(int number_subdivisions_input)
                                    :number_subdivisions(number_subdivisions_input)
    {
        backdoor_id=5;
        final_time=3;
        restart_step_number=0;
        output_directory="Simple_Ball/output";
        check_initial_mesh_for_self_intersection=false;
    }

    ~SIMPLE_BALL()
    {}

//#####################################################################
// Function Initialize_Cloth_State
//#####################################################################
void Initialize_Cloth_State(TRIANGLE_MESH*& triangle_mesh,PARTICLE_3D*& particles,TRIANGULATED_SURFACE*& triangulated_surface,
                                          DEFORMABLE_TRIANGULATED_SURFACE*& cloth)
{
    // initialize triangulated surface
    triangle_mesh=new TRIANGLE_MESH();
    particles=new PARTICLE_3D();
    triangulated_surface=new TRIANGULATED_SURFACE(*triangle_mesh,*particles);    
    cloth=new DEFORMABLE_TRIANGULATED_SURFACE(*triangulated_surface);
    if(restart_step_number) Read_Deformable_Triangulated_Surface(*cloth,restart_step_number);
    else{
        ARRAY<int> base_triangle(3,1,1);
        base_triangle.Set(1,1,2,3);
        triangle_mesh->Initialize_Triangle_Mesh(base_triangle);
        particles->Update_Position_And_Velocity();
        for(int k=0;k<triangle_mesh->number_nodes;k++) particles->array_collection->Add_Element();
        particles->X(1)=VECTOR_3D(0,.71,0);
        particles->X(2)=VECTOR_3D(2,.71,0);
        particles->X(3)=VECTOR_3D(1,.71,sqrt(3));
        for(int l=0;l<number_subdivisions;l++) triangulated_surface->Linearly_Subdivide();
        particles->Store_Mass();
        double mass_node=1./particles->number;
        copy(mass_node,particles->mass);}
}
//#####################################################################
// Function Initialize_Cloth_Dynamics
//#####################################################################
void Initialize_Cloth_Dynamics(DEFORMABLE_TRIANGULATED_SURFACE*& cloth)
{
    cloth->body_forces.Set_Gravity();
    cloth->body_forces.Set_Dynamic_Ether_Viscosity(.01);
    cloth->Set_CFL_Number(.5);
    cloth->Output_Artificial_Damping_Results();
    
    // springs
    cloth->triangulated_surface.triangle_mesh.Initialize_Segment_Mesh();
    UNIFORM_LINEAR_SPRINGS *springs=new UNIFORM_LINEAR_SPRINGS(*cloth->triangulated_surface.triangle_mesh.segment_mesh,cloth->triangulated_surface.particles);
    springs->Set_Stiffness(1);
    springs->Set_Restlength(2*pow(.5,number_subdivisions));
    springs->Set_Overdamping_Fraction(4);
     springs->artificially_damp_max_strain_per_time_step=false;
    springs->artificially_damp_strain=false;
    cloth->mesh_based_forces.Append(springs);

    // triangle_bending_springs
    TRIANGLE_BENDING_SPRINGS* bend=new TRIANGLE_BENDING_SPRINGS(cloth->triangulated_surface.particles);
    bend->Set_Quadruples_From_Triangle_Mesh(cloth->triangulated_surface.triangle_mesh);
    bend->Set_Stiffness(5);
    bend->Set_Damping(.0005);
    cloth->mesh_based_forces.Append(bend);

    // set up repulsion springs for collisions !!!!!!
    repulsion_springs_initialized=true;
    collisions_repulsion_spring_youngs_modulus=1;
    collisions_repulsion_spring_restlength=springs->restlength;
    collisions_repulsion_spring_mass=max(cloth->triangulated_surface.particles.mass);
}
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies(RIGID_BODY_LIST_3D<T>& solids_parameters.rigid_body_parameters.list)
{
    int index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/ground");
    index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/sphere",(T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->position=VECTOR_3D(1,.4,sqrt(3)/3);
    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Update_Object_Positions
//#####################################################################
void Update_Object_Positions(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const double time)
{
}
//#####################################################################
// Function Update_Object_Velocities
//#####################################################################
void Update_Object_Velocities(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const double time)
{
}
//#####################################################################
};
}
#endif

