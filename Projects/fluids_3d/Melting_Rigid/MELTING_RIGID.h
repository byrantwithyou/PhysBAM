//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_RIGID
//#####################################################################
#ifndef __MELTING_RIGID__
#define __MELTING_RIGID__

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include "../../solids/Embedded_Sphere/ANALYTIC_SPHERE.h"
#include "../../solids/Embedded_Torus/ANALYTIC_TORUS.h"
#include "../WATER_MELTING_EXAMPLE_3D.h"
#include <Heat_Flows/HEAT_3D.h>
#include <Level_Sets/EXTRAPOLATION_3D.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_3D.h>
#include <Rigid_Bodies/RIGID_BODY_LIST.h>
namespace PhysBAM{

template<class T,class RW=T>
class MELTING_RIGID:public WATER_MELTING_EXAMPLE_3D<T,RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::use_collision_aware_velocity_extrapolation;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::thin_shells_semi_lagrangian_phi;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::thin_shells_semi_lagrangian_velocity;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::use_collision_aware_signed_distance;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::clamp_phi_with_collision_bodies;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_substeps;
    using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::melting_parameters;using WATER_MELTING_EXAMPLE_3D<T,RW>::Construct_Levelsets_For_Objects;using WATER_MELTING_EXAMPLE_3D<T,RW>::maximum_velocity;
    using WATER_MELTING_EXAMPLE_3D<T,RW>::phi_objects;
    using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::initialized;using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::write_frame_title;
    using WATER_MELTING_EXAMPLE_3D<T,RW>::use_optimizations_for_constant_melting_speed;

    CYLINDER<T> cylinder_source;
    bool use_source;
    MATRIX<T,4> world_to_source;
    VECTOR<T,3> source_velocity;
    T source_start_time;

    // fluid parameters
    ANALYTIC_SPHERE<T> sphere;
    T initial_water_level;
    bool delete_positive_particles_crossing_bodies;

    // solid parameters
    int number_of_objects;
    ARRAY<VECTOR<T,3> > initial_position;
    ARRAY<QUATERNION<T> > initial_orientation;
    ARRAY<VECTOR<T,3> > initial_velocity;
    ARRAY<VECTOR<T,3> > initial_angular_velocity;

    LEVELSET_3D<GRID<TV> >* icecube_levelset;
    VECTOR<T,3> levelset_center;

    ARRAY<T,VECTOR<int,3> > pressure;
    ARRAY<T,VECTOR<int,3> > pressure_old;
    T ground_pressure_value;

    T melting_speed;
    T start_melting_time;

    bool use_two_way_coupling;

    T scaling;
    int m,n,mn;

    int example_number;

    TEMPERATURE_CONTAINER<T> temperature_container;
    bool use_temperature;
    T melting_temperature;
    T melting_gradient;

    bool rigid_body_preprocess;

    MELTING_RIGID(PARSE_ARGS& parse_args,int example_number_input)
        :use_two_way_coupling(true),scaling(.7),m(11),n(11),mn(11),temperature_container(fluids_parameters.grid_container),use_temperature(false)
    {
        example_number=example_number_input;
        std::cout<<"running MELTING_RIGID with example number "<<example_number<<std::endl;

        if(parse_args.Is_Value_Set("-use_temperature")) use_temperature=true;

        rigid_body_preprocess=false;
        if(parse_args.Is_Value_Set("-rigid_body_preprocess")) rigid_body_preprocess=true;

        if(example_number==3)use_optimizations_for_constant_melting_speed=true;

        // Two way coupling
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=true;
        solids_parameters.deformable_body_parameters.write=true;
        fluids_parameters.levelset_substeps=0.95;

        output_directory="Melting_Rigid/output_gen2";
        if(example_number==2)
            output_directory="Melting_Rigid/output_gen2-test";
        first_frame=0;last_frame=2000;frame_rate=48;
        restart=false;restart_frame=0;write_substeps=false;
        write_frame_title=true;
        //fluids_parameters.write_thin_shells_advection_cache=true;

        // Fluids parameters
        fluids_parameters.simulate=!rigid_body_preprocess;
        if(!fluids_parameters.simulate)
            output_directory="Melting_Rigid/rigid_body_preprocess";
        fluids_parameters.grid.Initialize(resolution,resolution,resolution,-1,1,0,2,-1,1);
        if(example_number==2)fluids_parameters.grid.Initialize(31,31,31,-1,1,0,2,-1,1);
        if(example_number==3)fluids_parameters.grid.Initialize(101,201,101,-1,1,0,4,-1,1);
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.reseeding_frame_rate=10;
        fluids_parameters.bias_towards_negative_particles=true;fluids_parameters.number_particles_per_cell=8;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=true;

        initial_water_level=-1;

        fluids_parameters.cfl=2.9;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;fluids_parameters.variable_viscosity=false;fluids_parameters.second_order_pressure=false;
        fluids_parameters.incompressible_iterations=100;

        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;
        fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.write_debug_data=false;fluids_parameters.restart_data_write_rate=1;
/*        fluids_parameters.write_levelset=false;fluids_parameters.write_velocity=false;
        fluids_parameters.write_particles=false;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=false;fluids_parameters.restart_data_write_rate=10000;*/

        fluids_parameters.store_particle_ids=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;

        // Solids parameters
        solids_parameters.gravity=9.8;
        solids_parameters.cfl=(T)1;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        verbose_dt=true;
        //solids_pre_roll=30;

        // Solids for melting parameters
        solids_parameters.collide_with_interior=false;
        solids_parameters.perform_self_collision=false;
        solids_parameters.synchronize_multiple_objects=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;

        // Melting parameters
        melting_parameters.maximum_depth=3;
        melting_parameters.refine_near_interface=true;
        melting_parameters.refine_for_high_deformation=true;

        melting_speed=10;
        number_of_objects=10;start_melting_time=2;ground_pressure_value=10;
        if(example_number==2)number_of_objects=1;
        if(example_number==3){number_of_objects=10;melting_speed=.03;}

        if(parse_args.Is_Value_Set("-number_of_objects")) number_of_objects=parse_args.Get_Integer_Value("-number_of_objects");
        if(parse_args.Is_Value_Set("-melting_speed")) melting_speed=(T)parse_args.Get_Double_Value("-melting_speed");
        if(parse_args.Is_Value_Set("-start_melting_time")) start_melting_time=(T)parse_args.Get_Double_Value("-start_melting_time");
        if(parse_args.Is_Value_Set("-ground_pressure_value")) ground_pressure_value=(T)parse_args.Get_Double_Value("-ground_pressure_value");

        initial_position.Resize(number_of_objects);
        initial_orientation.Resize(number_of_objects);
        initial_velocity.Resize(number_of_objects);
        initial_angular_velocity.Resize(number_of_objects);
//        basic_forces.Resize(number_of_objects);

        RANDOM_NUMBERS random;random.Set_Seed(1234);// set the seed so that the example is reproducible
        for(int i=0;i<number_of_objects;i++){
            initial_position(i)=VECTOR<T,3>(random.Get_Uniform_Number((T)-.4,(T).4),2.5+1.5*i,random.Get_Uniform_Number((T)-.4,(T).4));
            //initial_position(i)=VECTOR<T,3>(0,2.5+1.5*i,0);
            do{initial_orientation(i).s=random.Get_Uniform_Number(0,1);initial_orientation(i).v=random.Get_Uniform_Vector(VECTOR<T,3>(-1,-1,-1),VECTOR<T,3>(1,1,1));}while(initial_orientation(i).Magnitude()>1);
            initial_orientation(i).Normalize();
            initial_angular_velocity(i)=random.Get_Uniform_Vector(VECTOR<T,3>(-5,-5,-5),VECTOR<T,3>(5,5,5));}
        fluids_parameters.Initialize_Domain_Boundary_Conditions(); // sets up the proper wall states

        if(example_number==2){
            initial_position(1)=VECTOR<T,3>(.65,0.5,0);
            initial_orientation(1)=QUATERNION<T>(.25*pi,VECTOR<T,3>(0,1,0));
            initial_angular_velocity(1)=VECTOR<T,3>(0,0,0);
            //m=1;n=1;mn=1;melting_parameters.maximum_depth=1;
        }

        if(!rigid_body_preprocess&&!restart&&example_number==1){
            // read in the rigid body parameters from the other sim
//            /n/lie/disk2/data/losasso/Melting_Rigid_Slow/Melting_Rigid/output2
            std::cout<<"Reading rigid objects"<<std::endl;
            RIGID_BODY_LIST_3D<T> rigid_body_list;
            std::string restart_directory="/n/lie/disk2/data/losasso/Melting_Rigid_Slow/Melting_Rigid/output2";
            if(parse_args.Is_Value_Set("-restart_directory")) {
                std::cout<<"using non standard directory for restart"<<std::endl;
                restart_directory=parse_args.Get_String_Value("-restart_directory");}
            rigid_body_list.template Read<RW>(restart_directory,98);
            std::cout<<"-- Added "<<rigid_body_list.rigid_bodies.m<<" objects"<<std::endl;
            for(int i=6;i<=rigid_body_list.rigid_bodies.m;i++){
                RIGID_BODY<TV>& rigid_body=*rigid_body_list.rigid_bodies(i);//*rigid_body_list(i);
                initial_position(i-5)=rigid_body.position;
                initial_orientation(i-5)=rigid_body.orientation;
                rigid_body.Update_Angular_Velocity();
                initial_angular_velocity(i-5)=rigid_body.angular_velocity;
                initial_velocity(i-5)=rigid_body.velocity;
                }}


        icecube_levelset=new LEVELSET_3D<GRID<TV> >(*(new GRID<TV>()),*(new ARRAY<T,VECTOR<int,3> >()));

/*        FILE_UTILITIES::Read_From_File<RW>("/n/bcc/data/losasso/PhysBAM-Full/Public_Data/Rigid_Bodies/icecube.phi",*icecube_levelset);
        ARRAY<VECTOR<T,3> ,VECTOR<int,3> > V(icecube_levelset->grid,3);
        icecube_levelset->Set_Curvature_Motion(.01);
        T cfl=icecube_levelset->CFL(V);
        printf("\nSmoothing icecube:  ");
        int number_of_iterations=5;
        for(int i=0;i<5;i++)
        {
            printf("%d%% ",(int)((T)i*((T)100/(T)number_of_iterations)));fflush(stdout);
            for(int j=0;j<20;j++){printf(".");fflush(stdout);icecube_levelset->Euler_Step(V,cfl,0);}
            icecube_levelset->Reinitialize();
        }
        FILE_UTILITIES::Write_To_File<RW>("icecube.phi",*icecube_levelset);*/

        FILE_UTILITIES::Read_From_File<RW>("icecube.phi",*icecube_levelset);

        use_source=false;


        cylinder_source=CYLINDER<T>(VECTOR<T,3>(.65,1.9,0),VECTOR<T,3>(.65,2,0),.3);
        world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(0,-1,0);


        if(parse_args.Is_Value_Set("-other_source")){
            cylinder_source=CYLINDER<T>(VECTOR<T,3>(.95,1.5,-.5),VECTOR<T,3>(1.05,1.5,-.5),.3);
            world_to_source=MATRIX<T,4>::Identity_Matrix();source_velocity=VECTOR<T,3>(-1,0,0);}

        source_start_time=.5;
        if(parse_args.Is_Value_Set("-source_start_time")) source_start_time=(T)parse_args.Get_Double_Value("-source_start_time");

        if(use_temperature){
            temperature_container.Initialize_Array(3);
            melting_temperature=800;
            melting_gradient=melting_speed/(temperature_container.hot_point-melting_temperature);
            if(parse_args.Is_Value_Set("-melting_temperature")) melting_temperature=(T)parse_args.Get_Double_Value("-melting_temperature");
            temperature_container.Set_Velocity_3D(&fluids_parameters.incompressible.V);
        }
    }

    ~MELTING_RIGID()
    {}



//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Update_Fluid_Parameters(dt,time);
    if(time>=source_start_time) use_source=true;
}
//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    for(int i=0;i<number_of_objects;i++){
        Add_Melting_Object(melting_parameters.RIGID,0);}

    T domain_scale=1;
    T mu=(T).6;
    T epsilon=(T).2;
    if(example_number==3)
        epsilon=.8;

    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("/n/bcc/data/losasso/PhysBAM-Full/Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>("/n/bcc/data/losasso/PhysBAM-Full/Public_Data/Rigid_Bodies/ground",(T)domain_scale*.1);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
    solids_parameters.rigid_body_parameters.list(id)->position=domain_scale*VECTOR<T,3>(1,1,0);
    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(.5*pi,VECTOR<T,3>(0,0,1));
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>("/n/bcc/data/losasso/PhysBAM-Full/Public_Data/Rigid_Bodies/ground",(T)domain_scale*.1);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
    solids_parameters.rigid_body_parameters.list(id)->position=domain_scale*VECTOR<T,3>(-1,1,0);
    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(-.5*pi,VECTOR<T,3>(0,0,1));
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>("/n/bcc/data/losasso/PhysBAM-Full/Public_Data/Rigid_Bodies/ground",(T)domain_scale*.1);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
    solids_parameters.rigid_body_parameters.list(id)->position=domain_scale*VECTOR<T,3>(0,1,-1);
    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(.5*pi,VECTOR<T,3>(1,0,0));
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>("/n/bcc/data/losasso/PhysBAM-Full/Public_Data/Rigid_Bodies/ground",(T)domain_scale*.1);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction(mu);
    solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
    solids_parameters.rigid_body_parameters.list(id)->position=domain_scale*VECTOR<T,3>(0,1,1);
    solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(-.5*pi,VECTOR<T,3>(1,0,0));
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;

}
//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    ARRAY<RIGID_BODY<TV>*> extra_collision_bodies;
    Construct_Levelsets_For_Objects(time,extra_collision_bodies);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Grid(const int object,RED_GREEN_GRID_3D<T>& grid)
{
    printf("-------- %d ---------\n",object);
    GRID<TV>& levelset_grid=icecube_levelset->grid;
    levelset_center=levelset_grid.Domain().Center();
    T radius=scaling*.55*levelset_grid.Domain().Size().Max();
    grid.Initialize(GRID_3D<T>(4*m+1,4*n+1,4*mn+1,RANGE<T>(levelset_center).Thickened(radius)),melting_parameters.maximum_depth);
    sphere.center=grid.uniform_grid.Domain().Center();
    sphere.radius=.2;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    RED_GREEN_GRID_3D<T>& grid=melting_parameters.levelsets(object)->grid;
    ARRAY<VECTOR<T,3> >& node_locations=grid.Node_Locations();
    for(int p=0;p<phi.m;p++){
        phi(p)=icecube_levelset->Extended_Phi(((T)1/scaling)*(node_locations(p)-levelset_center)+levelset_center);}
/*
    VECTOR<T,3> center=grid.uniform_grid.Domain().Center();
    T size=.35;
    for(int p=0;p<phi.m;p++){
        VECTOR<T,3> distance=node_locations(p)-center;
        phi(p)=max(fabs(distance.x),fabs(distance.y),fabs(distance.z))-size;}*/
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR<T,3> >& V)
{
    //RED_GREEN_GRID_3D<T>& grid=melting_parameters.levelsets(1)->grid;
    //ARRAY<VECTOR<T,3> >& node_locations=grid.Node_Locations();

//    melting_parameters.use_constant_melting_speed=true;melting_parameters.constant_melting_speed=melting_speed;
    for(int p=0;p<V.m;p++) V(p)=VECTOR<T,3>(0,0,0);
//    for(int p=0;p<V.m;p++) V(p)=VECTOR<T,3>(node_locations(p).y-sphere.center.y,0,0);
//    for(int p=0;p<V.m;p++) V(p)=-node_locations(p)+circle.center;
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Melting_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(example_number==1||example_number==2){
        ARRAY<T,VECTOR<int,3> > smoothed_temperature=temperature_container.temperature_3d;
        HEAT_3D<T>::Smooth(fluids_parameters.grid,smoothed_temperature,3);
        //{ARRAY<T,VECTOR<int,3> > phi_object_temp(phi_objects);phi_object_temp*=(T)-1; // because extrapolation will extrapolate into positive phi
        //EXTRAPOLATION_3D<T,T> extrapolate(fluids_parameters.grid,phi_object_temp,smoothed_temperature);extrapolate.Set_Band_Width(15);extrapolate.Extrapolate();}

        T maximum_vel=0;
        for(int object=0;object<melting_parameters.levelsets.m;object++){
            LEVELSET_TETRAHEDRALIZED_VOLUME<T>& levelset=*melting_parameters.levelsets(object);
            RED_GREEN_GRID_3D<T>& grid=levelset.grid;
            int index=melting_parameters.body_index(object);
            if(!index)continue;
            RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(index);
            FRAME<T> frame=rigid_body.Frame()*melting_parameters.rigid_body_grid_frames(object).Inverse();
            // map data to the cell based indices
            ARRAY<VECTOR<T,3> >& node_locations=grid.Node_Locations();
            if(use_temperature)
            {
                LINEAR_INTERPOLATION<T,T> interpolation;
                for(int i=0;i<node_locations.m;i++){
                    T temperature=interpolation.Clamped_To_Array(fluids_parameters.grid,smoothed_temperature,frame*node_locations(i));
//                    if((frame*node_locations(i)).x>.6)melting_parameters.levelsets(object)->phi(i)+=dt;
                    if(temperature>melting_temperature)
                    {
                        T V=melting_gradient*(temperature-melting_temperature);
                        maximum_vel=max(maximum_velocity,V);
                        melting_parameters.levelsets(object)->phi(i)+=V*dt;
                    }
                }
            }
            else{
                maximum_vel=melting_speed;
                for(int i=0;i<node_locations.m;i++)if(fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.Phi(frame*node_locations(i))<0){
                    melting_parameters.levelsets(object)->phi(i)+=melting_speed*dt;}}}
        maximum_velocity=maximum_vel;
    }
    /*if(example_number==2){
        for(int object=0;object<melting_parameters.levelsets.m;object++){
            LEVELSET_TETRAHEDRALIZED_VOLUME<T>& levelset=*melting_parameters.levelsets(object);
            RED_GREEN_GRID_3D<T>& grid=levelset.grid;
            for(int i=0;i<grid.number_of_nodes;i++){
            melting_parameters.levelsets(object)->phi(i)+=.2*dt;}}}*/
    if(example_number==3)
        if(time>start_melting_time){
            melting_parameters.use_constant_melting_speed=true;melting_parameters.constant_melting_speed=min((time-2)*(T).25,melting_speed);}

    WATER_MELTING_EXAMPLE_3D<T,RW>::Melting_Substep(dt,time);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Particle_Positions_And_Velocities(const int object)
{
    RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(melting_parameters.body_index(object));
    //rigid_body.is_static=true;
    rigid_body.position=initial_position(object);
    rigid_body.velocity=initial_velocity(object);
    rigid_body.orientation=initial_orientation(object);
    rigid_body.angular_velocity=initial_angular_velocity(object);
    rigid_body.Update_Angular_Momentum();
}
//#####################################################################
// Function Initialize_Forces
//#####################################################################
void Initialize_Forces()
{
    PHYSBAM_FATAL_ERROR("BASIC FORCES NO LONGER EXIST");
#if 0
    for(int object=0;object<melting_parameters.body_index.m;object++){
        int index=melting_parameters.body_index(object);if(!index)continue;
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(index);
        std::cout<<rigid_body.position<<std::endl;
        rigid_body.Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

        if(!rigid_body_preprocess&&use_two_way_coupling){
            basic_forces(object)=new RIGID_BODY_BASIC_FORCES_3D<T>(rigid_body);
            basic_forces(object)->V_grid=fluids_parameters.incompressible.projection.p_grid;
//            basic_forces->Set_Wind_Pressure(pressure,3); // 50%
            basic_forces(object)->Set_Wind_Pressure(pressure,2); // floats pretty high
            basic_forces(object)->Use_Wind_Drag();
            rigid_body.body_forces.Append(basic_forces(object));}
        /*if(use_two_way_coupling){
            static int count=0;count++;
            Add_To_Fluid_Simulation(rigid_body,true,false,.1);}*/
    }
#endif
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    temperature_container.Initialize_Array(3);
    // Not so good to set up a heaviside function here because then the interface will
    // be exactly between the two nodes which can lead to roundoff issues when setting dirichlet cells, etc.
    GRID<TV>& grid=fluids_parameters.grid;
    fluids_parameters.particle_levelset_evolution.particle_levelset.Set_Minimum_Particle_Radius((T).2*grid.max_dx_dy_dz);
    fluids_parameters.particle_levelset_evolution.particle_levelset.Set_Maximum_Particle_Radius((T).7*grid.max_dx_dy_dz);
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++)
        fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=grid.y(j)-grid.ymin-initial_water_level;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
template<class SOURCE> void Adjust_Phi_With_Sources(const SOURCE& source,const T time)
{
    if(!use_source) return;
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++){
        VECTOR<T,3> source_X=world_to_source*grid.X(i,j,ij);
        if(source.Lazy_Inside(source_X)) fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=min(fluids_parameters.particle_levelset_evolution.phi(i,j,ij),source.Signed_Distance(source_X));}
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    if(!use_source) return;
    Adjust_Phi_With_Sources(cylinder_source,time);
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
template<class SOURCE> void Get_Source_Reseed_Mask(const SOURCE& source,ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,3> >(grid);
    T padding=3*grid.max_dx_dy_dz;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) if(!source.Outside(world_to_source*grid.X(i,j,ij),padding)) (*cell_centered_mask)(i,j,ij)=true;
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time)
{
    if(!use_source) return;
    Get_Source_Reseed_Mask(cylinder_source,cell_centered_mask,time);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class SOURCE> void Get_Source_Velocities(const SOURCE& source,const T time)
{
    GRID<TV> &u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid,&w_grid=fluids_parameters.w_grid;
    PROJECTION_3D<T>& projection=fluids_parameters.incompressible.projection;
    for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) for(int ij=0;ij<u_grid.mn;ij++) if(source.Lazy_Inside(world_to_source*u_grid.X(i,j,ij))){projection.elliptic_solver->psi_N_u(i,j,ij)=true;projection.u(i,j,ij)=source_velocity.x;}
    for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) for(int ij=0;ij<v_grid.mn;ij++) if(source.Lazy_Inside(world_to_source*v_grid.X(i,j,ij))){projection.elliptic_solver->psi_N_v(i,j,ij)=true;projection.v(i,j,ij)=source_velocity.y;}
    for(int i=0;i<w_grid.m;i++) for(int j=0;j<w_grid.n;j++) for(int ij=0;ij<w_grid.mn;ij++) if(source.Lazy_Inside(world_to_source*w_grid.X(i,j,ij))){projection.elliptic_solver->psi_N_w(i,j,ij)=true;projection.w(i,j,ij)=source_velocity.z;}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    if(!use_source) return;
    Get_Source_Velocities(cylinder_source,time);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Update_Melting_Substep_Parameters(const T dt,const T time)
{
    if(use_temperature){
        GRID<TV>& grid=fluids_parameters.grid;
        temperature_container.Euler_Step(dt,time);
        for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++){
            if(fluids_parameters.particle_levelset_evolution.phi(i,j,ij)<0){
                if(use_source){
                    if(cylinder_source.Lazy_Inside(world_to_source*grid.X(i,j,ij)))
                        temperature_container.temperature_3d(i,j,ij)=temperature_container.hot_point;}}
            /*else temperature_container.temperature_3d(i,j,ij)=temperature_container.ambient_temperature;*/}}

    for(int i=0;i<basic_forces.m;i++)if(basic_forces(i)){
        basic_forces(i)->V_grid=fluids_parameters.incompressible.projection.p_grid;
        //basic_forces(i)->wind_pressure_scaling_factor=0;
    }
    static bool first=true;
    pressure_old=pressure;
    pressure=fluids_parameters.incompressible.projection.p;

    GRID<TV>& grid=fluids_parameters.incompressible.projection.p_grid;
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)
        if(phi_objects(i,j,ij)+phi_objects(i+1,j,ij)+phi_objects(i,j+1,ij)+phi_objects(i+1,j+1,ij)+phi_objects(i,j,ij+1)+phi_objects(i+1,j,ij+1)+phi_objects(i,j+1,ij+1)+phi_objects(i+1,j+1,ij+1)<0)
            pressure(i,j,ij)=0;
    /*if(time>start_melting_time)*/
    for(int i=0;i<grid.m;i++)for(int j=0;j<=2;j++)for(int ij=0;ij<grid.mn;ij++)
        pressure(i,j,ij)=10;//min((T)ground_pressure_value,(time-2)*ground_pressure_value);
    for(int i=0;i<=2;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)
        pressure(i,j,ij)=5;//min((T)ground_pressure_value,(time-2)*ground_pressure_value);
    for(int i=grid.m-2;i<=grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)
        pressure(i,j,ij)=5;//min((T)ground_pressure_value,(time-2)*ground_pressure_value);
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<=2;ij++)
        pressure(i,j,ij)=5;//min((T)ground_pressure_value,(time-2)*ground_pressure_value);
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=grid.mn-2;ij<=grid.mn;ij++)
        pressure(i,j,ij)=5;//min((T)ground_pressure_value,(time-2)*ground_pressure_value);

    printf("\n\nMINIMUM PRESSURE IS: %f ---- MAXIMUM PRESSURE IS: %f\n",pressure.Min(),pressure.Max());

    HEAT_3D<T>::Smooth(grid,pressure,5);
    if(!first)ARRAY<T,VECTOR<int,3> >::copy(.5,pressure,.5,pressure_old,pressure);
    printf("\n\nMINIMUM PRESSURE IS: %f ---- MAXIMUM PRESSURE IS: %f\n",pressure.Min(),pressure.Max());
    first=false;
//    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)
//        pressure(i,j,ij)=1000*(4-grid.X(i,j,ij).y);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    WATER_MELTING_EXAMPLE_3D<T,RW>::Read_Output_Files_Solids(frame);
    if(!rigid_body_preprocess&&use_temperature){
        std::string prefix=output_directory+"/";
        FILE_UTILITIES::Read_From_File<RW>(prefix+STRING_UTILITIES::string_sprintf("temperature.%d",frame),temperature_container.temperature_3d);
    }
    fluids_parameters.particle_levelset_evolution.particle_levelset.Reseed_Particles();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    WATER_MELTING_EXAMPLE_3D<T,RW>::Write_Output_Files(frame);
    if(!rigid_body_preprocess&&use_temperature){
        std::string prefix=output_directory+"/";
        FILE_UTILITIES::Write_To_File<RW>(prefix+STRING_UTILITIES::string_sprintf("temperature.%d",frame),temperature_container.temperature_3d);}
}

//#####################################################################
};
}
#endif


