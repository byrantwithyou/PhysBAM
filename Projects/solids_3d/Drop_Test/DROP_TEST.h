//#####################################################################
// Copyright 2002, Robert Bridson, Ronald Fedkiw, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DROP_TEST
//##################################################################### 
//
// Arbitrary triangulated surface dropped on the ground
//
//#####################################################################
// Bridson - June 27, 2002
//#####################################################################
#ifndef __DROP_TEST__
#define __DROP_TEST__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <fstream>
#include "../CLOTH_EXAMPLE.h"
namespace PhysBAM{

class DROP_TEST:public CLOTH_EXAMPLE
{
public:
    TRIANGLE_MESH initial_mesh;
    PARTICLE_3D initial_particles;
    TRIANGULATED_SURFACE initial_surface;

    explicit DROP_TEST(const char *surface_file_name)
                                    :initial_surface(initial_mesh,initial_particles)
    {
        std::ifstream input(surface_file_name,std::ios::binary);
        initial_surface.Read(input);
        input.close();
        double min_y=initial_particles.X(1).y;
        int i;
        for(i=2; i<=initial_particles.array_collection->Size(); ++i) min_y=min(initial_particles.X(i).y,min_y);
        for(i=1; i<=initial_particles.array_collection->Size(); ++i) initial_particles.X(i).y+=1e-3-min_y;
        backdoor_id=4;
        final_time=10;
        restart_step_number=0;
        output_directory="Drop_Test/output";
    }

    ~DROP_TEST()
    {}

//#####################################################################
// Function Initialize_Cloth_State
//#####################################################################
void Initialize_Cloth_State(TRIANGLE_MESH*& triangle_mesh,PARTICLE_3D*& particles,TRIANGULATED_SURFACE*& triangulated_surface,
                                          DEFORMABLE_TRIANGULATED_SURFACE*& cloth)
{
    triangle_mesh=new TRIANGLE_MESH(initial_mesh);
    particles=new PARTICLE_3D(initial_particles);
    particles->Update_Position_And_Velocity();
    particles->Store_Mass();
    for(int t=1; t<=triangle_mesh->triangles.m; ++t){
        int i,j,k;triangle_mesh->triangles.Get(t,i,j,k);
        TRIANGLE_3D tri(particles->X(i),particles->X(j),particles->X(k));
        double triangle_mass=tri.Area();
        particles->mass(i)+=triangle_mass/3;
        particles->mass(j)+=triangle_mass/3;
        particles->mass(k)+=triangle_mass/3;
    }
    triangulated_surface=new TRIANGULATED_SURFACE(*triangle_mesh,*particles);
    cloth=new DEFORMABLE_TRIANGULATED_SURFACE(*triangulated_surface);
    
    if(restart_step_number) Read_Deformable_Triangulated_Surface(*cloth,restart_step_number);
}
//#####################################################################
// Function Initialize_Cloth_Dynamics
//#####################################################################
void Initialize_Cloth_Dynamics(DEFORMABLE_TRIANGULATED_SURFACE*& cloth)
{
    cloth->body_forces.Set_Gravity();
    cloth->Set_CFL_Number(.9);
    cloth->Output_Artificial_Damping_Results();
    
    // springs
    initial_mesh.Initialize_Segment_Mesh();
    LINEAR_SPRINGS *springs=new LINEAR_SPRINGS(*initial_mesh.segment_mesh,cloth->triangulated_surface.particles);
    springs->Set_Stiffness(1);
    if(restart_step_number){
        ARRAY<VECTOR_3D,VECTOR<int,1> > Xsave(cloth->triangulated_surface.particles.X);
        cloth->triangulated_surface.particles.X=initial_particles.X;
        springs->Set_Restlength_From_Particles();
        cloth->triangulated_surface.particles.X=Xsave;
    }else springs->Set_Restlength_From_Particles();
    springs->Set_Overdamping_Fraction(2); 
    springs->Set_Artificial_Max_Strain_Per_Time_Step(.1);
    springs->Set_Artificial_Min_And_Max_Strain(-.025,.15);
    cloth->mesh_based_forces.Append(springs);

    // triangle_bending_springs
    TRIANGLE_BENDING_SPRINGS* bend=new TRIANGLE_BENDING_SPRINGS(cloth->triangulated_surface.particles);
    bend->Set_Quadruples_From_Triangle_Mesh(initial_mesh);
    if(restart_step_number){
        ARRAY<VECTOR_3D,VECTOR<int,1> > Xsave(cloth->triangulated_surface.particles.X);
        cloth->triangulated_surface.particles.X=initial_particles.X;
        bend->Set_Constants_From_Particles(1,1e-2);
        cloth->triangulated_surface.particles.X=Xsave;
    }else bend->Set_Constants_From_Particles(1,1e-2);
    cloth->mesh_based_forces.Append(bend);

    // set up repulsion springs for collisions !!!!!!
    repulsion_springs_initialized=true;
    collisions_repulsion_spring_youngs_modulus=1;
    collisions_repulsion_spring_restlength=max(springs->restlength);
    collisions_repulsion_spring_mass=max(cloth->triangulated_surface.particles.mass);
}
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies(RIGID_BODY_LIST_3D<T>& solids_parameters.rigid_body_parameters.list)
{
    int index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/ground");
/*
    // ceiling for crushing
    rigid_bodies.Append(new RIGID_BODY<TV>);
    rigid_bodies_friction.Append(1);
    // triangulated surface
    index=triangulated_surface_list.Add_Triangulated_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.tri");
    input.open(filename,std::ios::in|std::ios::binary);
    triangulated_surface_list.triangulated_surface(index)->Read(input);input.close();
    rigid_bodies(rigid_bodies.m)->Initialize_Triangulated_Surface(*triangulated_surface_list.triangulated_surface(index));
    // implicit surface
    index=implicit_surface_list.Add_Levelset_Implicit_Surface();
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.phi");
    input.open(filename,std::ios::in|std::ios::binary);
    implicit_surface_list.implicit_surface(index)->Read(input);input.close();
    rigid_bodies(rigid_bodies.m)->Initialize_Implicit_Surface(*implicit_surface_list.implicit_surface(index));
    // rigid body
    sprintf(filename,"../../Public_Data/Rigid_Bodies/ground.rgd");
    input.open(filename,std::ios::in|std::ios::binary);
    rigid_bodies(rigid_bodies.m)->Read(input);input.close();
    // transform
    rigid_bodies(rigid_bodies.m)->position=VECTOR_3D(0,.95,0);
    rigid_bodies(rigid_bodies.m)->orientation=QUATERNION(pi,VECTOR_3D(1,0,0))
*/
    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Update_Object_Positions
//#####################################################################
void Update_Object_Positions(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const double time)
{
    /*
    rigid_bodies(2)->position=VECTOR_3D(0,max(.95-time/3,.05),0);
    rigid_bodies(2)->orientation=QUATERNION(pi,VECTOR_3D(1,0,0));
    */
}
//#####################################################################
// Function Update_Object_Velocities
//#####################################################################
void Update_Object_Velocities(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const double time)
{
    /*
    if(.95-time>.05)
        rigid_bodies(2)->velocity=VECTOR_3D(0,-1/3,0);
    else
        rigid_bodies(2)->velocity=VECTOR_3D(0,0,0);
    */
}
//#####################################################################
};
}
#endif

