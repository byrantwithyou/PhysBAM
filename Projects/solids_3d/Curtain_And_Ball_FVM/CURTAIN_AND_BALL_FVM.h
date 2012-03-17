//#####################################################################
// Copyright 2002, 2003, Robert Bridson, Ronald Fedkiw, Neil Molino, Zhaosheng Bao.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CURTAIN_AND_BALL_FVM
//##################################################################### 
// Hanging curtain with a moving ball
#ifndef __CURTAIN_AND_BALL_FVM__
#define __CURTAIN_AND_BALL_FVM__

#include "../CLOTH_EXAMPLE.h"
namespace PhysBAM{

template<class T>
class CURTAIN_AND_BALL_FVM:public CLOTH_EXAMPLE<T>
{
public:
    int number_side_panels;
    T aspect_ratio,side_length;

    CURTAIN_AND_BALL_FVM()
        :number_side_panels(40),aspect_ratio((T)1),side_length(1)
    {
        final_time=7;
        restart_step_number=0;
        cfl_number=(T)1;
        cg_tolerance=(T)1e-3;
        cg_iterations=1000;
        
        // Springs Stuff for Comparison
        use_masses_and_springs=false;use_altitude_springs=false;
        edge_spring_stiffness=(T)1e5;edge_spring_overdamping_fraction=2;
        use_shortest_altitude_spring_only=true;
        use_altitude_springs_compressed_by_threshold_only=true;altitude_spring_fraction_compression=(T).1;
        altitude_spring_stiffness=2*edge_spring_stiffness;altitude_spring_overdamping_fraction=2*edge_spring_overdamping_fraction;
        
        // Finite Volumes
        use_fvm=false;
        use_diagonalized_fvm=true;
        youngs_modulus=100;
        poissons_ratio=(T).3;
        Rayleigh_coefficient=(T).02;
        use_linear_elasticity=false;

        // Bending Forces
        bending_stiffness=(T).1;bending_damping=(T).02*bending_stiffness;
        use_bending_plasticity=true;
        bending_plastic_yield=1;
        bending_plastic_hardening=1;
        
        // Gravity
        output_directory="Curtain_And_Ball_FVM/output";
        check_initial_mesh_for_self_intersection=false;
        perform_self_collision=true;
        gravity=10;
    }
  
    ~CURTAIN_AND_BALL_FVM()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data(TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    
    particles.Update_Position_And_Velocity();particles.Store_Mass();

    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    triangle_mesh.Initialize_Herring_Bone_Mesh(m,n);
    for(int k=0;k<triangle_mesh.number_nodes;k++) particles.Add_Element();
    T mass_node=aspect_ratio*sqr(side_length)/(m*n);ARRAY<T,VECTOR<int,1> >::copy(mass_node,particles.mass);    
    T dx=aspect_ratio*side_length/(m-1),dy=side_length/(n-1);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){int node=i+m*(j-1);
        particles.X(node)=VECTOR_3D<T>((i-1)*dx,.5,(j-1)*dy);
        particles.V(node)=VECTOR_3D<T>(0,0,0);}
}
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies()
{
    int index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    index=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body("../../Public_Data/Rigid_Bodies/sphere",(T).25);
    rigid_bodiy_list.rigid_bodies(index)->position.z+=.5; // reset position
    rigid_bodiy_list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T)0);
    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Update_Rigid_Body_Positions
//#####################################################################
void Update_Rigid_Body_Positions(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const T time)
{
    if(time < 2) rigid_bodies(2)->position=VECTOR_3D<T>(0,0,.5);
    else if(time < 3.5) rigid_bodies(2)->position=VECTOR_3D<T>((time-2),(T).5*(time-2),.5);
    else if(time < 4) rigid_bodies(2)->position=VECTOR_3D<T>(1.5,(T)(.75-1.5*(time-3.5)),.5);
    else rigid_bodies(2)->position=VECTOR_3D<T>((T)(1.5-1.5*(time-4)),0,.5);
}
//#####################################################################
// Function Update_Rigid_Body_Velocities
//#####################################################################
void Update_Rigid_Body_Velocities(ARRAY<RIGID_BODY<TV>*>& rigid_bodies,const T time)
{
    if(time < 2) rigid_bodies(2)->velocity=VECTOR_3D<T>(0,0,0);
    else if(time < 3.5) rigid_bodies(2)->velocity=VECTOR_3D<T>(1,.5,0);
    else if(time < 4) rigid_bodies(2)->velocity=VECTOR_3D<T>(0,-1.5,0);
    else rigid_bodies(2)->velocity=VECTOR_3D<T>(-1.5,0,0);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
// for external forces and velocities
void Set_External_Velocities(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
      int i,j;
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    i=1;j=1;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
    i=1;j=n;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
    //i=m;j=n;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
    //i=m;j=1;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
// for external forces and velocities
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    int i,j;
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    i=1;j=1;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
    i=1;j=n;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
    //i=m;j=n;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
    //i=m;j=1;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
};
}
#endif
