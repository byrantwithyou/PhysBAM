//#####################################################################
// Copyright 2002-2004, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CURTAIN_AND_BALL
//##################################################################### 
// Hanging curtain with a moving ball
#ifndef __CURTAIN_AND_BALL__
#define __CURTAIN_AND_BALL__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Forces_And_Torques/BODY_FORCES_3D.h>
namespace PhysBAM{

template<class T,class RW>
class CURTAIN_AND_BALL:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;
    using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;

    int number_side_panels;
    T aspect_ratio,side_length;
    bool collide_against_sphere_triangulated_surface,collide_against_sphere_particles;
    RIGID_BODY_LIST_3D<T> rigid_bodies_to_collide_against;
    PARTICLES<T,VECTOR_3D<T> > sphere_particles;
    TRIANGULATED_SURFACE<T>* sphere_surface;

    CURTAIN_AND_BALL()
        :BASE(0,fluids_parameters.NONE),number_side_panels(40),aspect_ratio((T)1.7),side_length(1),collide_against_sphere_triangulated_surface(false),
         collide_against_sphere_particles(false),sphere_surface(0)
    {
        //allow_intersections=true;allow_intersections_tolerance=(T)1e-3;
        last_frame=7*24;
        restart=false;restart_frame=0;
        solids_parameters.cfl=(T)5.9;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        output_directory="Curtain_And_Ball/output";
        //min_collision_loops=max_collision_loops=1;
    }

    ~CURTAIN_AND_BALL()
    {delete sphere_surface;}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).triangulated_surface;
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;

    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    triangle_mesh.Initialize_Herring_Bone_Mesh(m,n);
    for(int k=1;k<=triangle_mesh.number_nodes;k++) particles.array_collection->Add_Element();
    T mass_node=aspect_ratio*sqr(side_length)/(m*n);ARRAY<T>::copy(mass_node,particles.mass.array); 
    T dx=aspect_ratio*side_length/(m-1),dy=side_length/(n-1);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){int node=i+m*(j-1);
        particles.X(node)=VECTOR_3D<T>((i-1)*dx,.5,(j-1)*dy);
        particles.V(node)=VECTOR_3D<T>(0,0,0);}

    //rigid bodies
    int rigid_index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index)->is_static=true;
    rigid_bodies_to_collide_against.rigid_bodies.Append(solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index));
    rigid_index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>("../../Public_Data/Rigid_Bodies/sphere",(T).25);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index)->frame.t.z+=.5; // reset position
    solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index)->Set_Coefficient_Of_Friction((T)0);

    if(collide_against_sphere_triangulated_surface||collide_against_sphere_particles){
        sphere_particles.Update_Velocity();sphere_particles.Store_Mass();int sphere_index=2;
        for(int p=1;p<=solids_parameters.rigid_body_parameters.list.rigid_bodies(sphere_index)->triangulated_surface->particles.X.array.m;p++){
            sphere_particles.array_collection->Add_Element();sphere_particles.X(p)=solids_parameters.rigid_body_parameters.list.rigid_bodies(sphere_index)->triangulated_surface->particles.X(p);
            sphere_particles.V(p)=solids_parameters.rigid_body_parameters.list.rigid_bodies(sphere_index)->velocity;
            sphere_particles.mass(p)=(T)1e6;}
        if(collide_against_sphere_triangulated_surface){
            sphere_surface=new TRIANGULATED_SURFACE<T>(solids_parameters.rigid_body_parameters.list.rigid_bodies(sphere_index)->triangulated_surface->triangle_mesh,sphere_particles);
            solid_body_collection.deformable_body_collection.triangle_collisions.geometry.Add_Triangulated_Surface(*sphere_surface);}
        else if(collide_against_sphere_particles) solid_body_collection.deformable_body_collection.triangle_collisions.geometry.Add_Free_Particles(sphere_particles);
        solids_parameters.collision_body_list.Add_Bodies(rigid_bodies_to_collide_against);}
    else solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Deformable_Objects
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();

    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Body_Forces<T>(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface));
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Edge_Springs<T>(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface,2/(1+sqrt((T)2)),2));
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Altitude_Springs<T>(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface,2*4/(1+sqrt((T)2)),4));
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Bending_Elements(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface));
    //solids_parameters.deformable_body_parameters.list(1).Add_Bending_Springs(solids_parameters.deformable_body_parameters.list(1).triangulated_surface->triangle_mesh);
}
//#####################################################################
// Function Update_Collision_Body_Positions_And_Velocities
//#####################################################################
void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE
{
    if(time < 2){solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->frame.t=VECTOR_3D<T>(0,0,.5);solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->velocity=VECTOR_3D<T>(0,0,0);}
    else if(time < 3.5){solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->frame.t=VECTOR_3D<T>((time-2),(T).5*(time-2),.5);solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->velocity=VECTOR_3D<T>(1,.5,0);}
    else if(time < 4){solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->frame.t=VECTOR_3D<T>(1.5,(T)(.75-1.5*(time-3.5)),.5);solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->velocity=VECTOR_3D<T>(0,-1.5,0);}
    else{solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->frame.t=VECTOR_3D<T>((T)(1.5-1.5*(time-4)),0,.5);solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->velocity=VECTOR_3D<T>(-1.5,0,0);}
    if(collide_against_sphere_triangulated_surface||collide_against_sphere_particles){
        FRAME_3D<T> frame=solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->frame;
        for(int p=1;p<=sphere_particles.X.array.m;p++){
            sphere_particles.X(p)=frame*solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->triangulated_surface->particles.X(p);
            sphere_particles.V(p)=solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->velocity;}}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
// for external forces and velocities
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time) PHYSBAM_OVERRIDE
{
    switch(id_number){
    case 1:
        {int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=1;j=1;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);i=1;j=n;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
// for external forces and velocities
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time) PHYSBAM_OVERRIDE
{
    switch(id_number){
    case 1:
        {int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=1;j=1;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);i=1;j=n;V(i+m*(j-1))=VECTOR_3D<T>(0,0,0);}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
};
}
#endif

