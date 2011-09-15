//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CLOTH_INCLINE__
#define __CLOTH_INCLINE__

#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <Forces_And_Torques/BODY_FORCES_3D.h>
namespace PhysBAM{

template<class T,class RW>
class CLOTH_INCLINE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::frame_rate;
    using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;using BASE::initial_time;

    int number_side_panels;
    T aspect_ratio,side_length;
    T initial_cloth_height;
    FRAME_3D<T> incline_frame;
    T ground_friction_coefficient;
    enum {LEVELSET,COARSE_TRIANGULATED_SURFACE,DENSE_TRIANGULATED_SURFACE} collision_object_type;
    PARTICLES<T,VECTOR_3D<T> > plane_particles; // transformed plane particles
    TRIANGULATED_SURFACE<T>* plane_surface;
    bool use_repulsions_only;
    
    CLOTH_INCLINE()
        :BASE(0,fluids_parameters.NONE),aspect_ratio((T)1.5),side_length(1),plane_surface(0),use_repulsions_only(true)
    {
        last_frame=(int)(3*frame_rate); // 7 seconds
        restart=false;restart_frame=0;
        solids_parameters.cfl=(T)5.9;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;

        PARAMETER_LIST parameters;
        parameters.Begin_Parse("sim.params");
        number_side_panels=parameters.Get_Parameter("number_side_panels",(int)40);
        output_directory=parameters.Get_Parameter("output_directory",std::string("Cloth_Incline/output"));
        T incline_angle_degrees=parameters.Get_Parameter("incline_angle_degrees",(T)30);
        incline_frame=FRAME_3D<T>(VECTOR_3D<T>(),QUATERNION<T>(pi/(T)180*incline_angle_degrees,VECTOR_3D<T>(1,0,0)));
        solids_parameters.self_collision_friction_coefficient=ground_friction_coefficient=parameters.Get_Parameter("ground_friction_coefficient",(T).4);
        initial_cloth_height=parameters.Get_Parameter("initial_cloth_height",(T).02);
        solids_parameters.collisions_repulsion_thickness=parameters.Get_Parameter("collisions_repulsion_thickness",(T)1e-3);
        solids_parameters.turn_off_all_collisions=parameters.Get_Parameter("turn_off_all_collisions",(bool)true);
        std::string collision_type=parameters.Get_Parameter("collision_object_type",std::string("DENSE_TRIANGULATED_SURFACE"));
        if(collision_type=="LEVELSET") collision_object_type=LEVELSET;
        else if(collision_type=="DENSE_TRIANGULATED_SURFACE") collision_object_type=DENSE_TRIANGULATED_SURFACE;
        else if(collision_type=="COARSE_TRIANGULATED_SURFACE") collision_object_type=COARSE_TRIANGULATED_SURFACE;
        parameters.End_Parse();
        parameters.Write(output_directory+"/sim.params");
    }

    //#####################################################################
    // Function Get_Initial_Data
    //#####################################################################
    void Get_Initial_Data()
    {int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).triangulated_surface;
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;
    
    // initialize cloth
    GRID<TV> cloth_grid(number_side_panels+1,(int)(aspect_ratio*number_side_panels)+1,0,1,0,aspect_ratio);
    triangle_mesh.Initialize_Herring_Bone_Mesh(cloth_grid.m,cloth_grid.n);
    particles.array_collection->Add_Elements(triangle_mesh.number_nodes);
    for(int i=1;i<=cloth_grid.m;i++) for(int j=1;j<=cloth_grid.n;j++){
        int node=i+cloth_grid.m*(j-1);
        particles.X(node)=VECTOR_3D<T>(cloth_grid.x(i),initial_cloth_height,cloth_grid.y(j));
        particles.V(node)=VECTOR_3D<T>();}
    T mass_node=aspect_ratio*sqr(side_length)/(cloth_grid.m*cloth_grid.n);ARRAY<T>::copy(mass_node,particles.mass.array); 
    
    // read collision object
    int rigid_index=0;
    if(collision_object_type==LEVELSET) rigid_index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>(data_directory+"/Rigid_Bodies/ground");
    else if(collision_object_type==COARSE_TRIANGULATED_SURFACE) rigid_index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Triangulated_Surfaces/small_plane_1k",1,true,false,false,false);
    else if(collision_object_type==DENSE_TRIANGULATED_SURFACE) rigid_index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Triangulated_Surfaces/small_plane_20k",1,true,false,false,false);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index)->Set_Coefficient_Of_Friction(ground_friction_coefficient);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index)->is_static=true;
    solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index)->frame=incline_frame;

    // set collision method
    if(collision_object_type==LEVELSET){
        std::cout<<"Using Level Set Collisions"<<std::endl;
        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);}
    else if(collision_object_type==COARSE_TRIANGULATED_SURFACE || collision_object_type==DENSE_TRIANGULATED_SURFACE){
        std::cout<<"Using triangle collisions"<<std::endl;
        plane_particles.Update_Velocity();plane_particles.Store_Mass();
        plane_particles.array_collection->Add_Elements(solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index)->triangulated_surface->particles.array_collection->Size());
        for(int p=1;p<=plane_particles.array_collection->Size();p++) plane_particles.mass(p)=(T)1e6;
        plane_surface=new TRIANGULATED_SURFACE<T>(solids_parameters.rigid_body_parameters.list.rigid_bodies(rigid_index)->triangulated_surface->triangle_mesh,plane_particles);
        Update_Collision_Body_Positions_And_Velocities(initial_time);
        solids_parameters.extra_collision_surfaces.Append(plane_surface);}}
    //#####################################################################
    // Function Initialize_Bodies
    //#####################################################################
    void Initialize_Bodies() PHYSBAM_OVERRIDE
    {Get_Initial_Data();
    // setup springs
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Body_Forces<T>(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface));
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Edge_Springs<T>(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface,2/(1+sqrt((T)2)),2));
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Altitude_Springs<T>(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface,2*4/(1+sqrt((T)2)),4));
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Bending_Elements(*solids_parameters.deformable_body_parameters.list(1).triangulated_surface));}
    //#####################################################################
    // Function Update_Collision_Body_Positions_And_Velocities
    //#####################################################################
    void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE
    {if(collision_object_type==COARSE_TRIANGULATED_SURFACE || collision_object_type==DENSE_TRIANGULATED_SURFACE){
        const RIGID_BODY<TV>& body=*solids_parameters.rigid_body_parameters.list.rigid_bodies(1);
        for(int p=1;p<=plane_particles.array_collection->Size();p++){
            plane_particles.X(p)=body.World_Space_Point(body.triangulated_surface->particles.X(p));
            plane_particles.V(p)=body.Pointwise_Object_Velocity(plane_particles.X(p));}}}

//#####################################################################
};
}
#endif
