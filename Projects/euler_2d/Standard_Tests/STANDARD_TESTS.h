//#####################################################################
// Copyright 2009 Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include "math.h"
#include <fstream>
#include <iostream>

#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Rigids/Rigid_Bodies/MASS_PROPERTIES.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>

namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
public:
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    typedef VECTOR<T,2*TV::m> T_FACE_VECTOR;typedef VECTOR<TV,2*TV::m> TV_FACE_VECTOR;

    using BASE::last_frame;using BASE::frame_rate;using BASE::viewer_dir;using BASE::fluids_parameters;using BASE::solids_parameters;
    using BASE::solids_fluids_parameters;using BASE::Add_To_Fluid_Simulation;
    using BASE::stream_type;using BASE::data_directory;
    using BASE::solid_body_collection;using BASE::test_number;using BASE::resolution;
    using BASE::user_last_frame;

    SOLIDS_STANDARD_TESTS<TV> tests;
    SEGMENTED_CURVE_2D<T>* boundary;
    ARRAY<int> bound_particles;
    std::ofstream gnuplot_file_stream;
    T spring_factor;
    TV original_position;

    TV_DIMENSION state_left,state_right; // // (density,velocity_x,velocity_y,pressure)
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    int eno_scheme;
    int eno_order;
    int rk_order;
    T cfl_number;
    bool timesplit;
    bool exact;
    bool use_slip;

    /***************
    example explanation:
    1. Cylinder lift-off.
    2. Inflating bladder.
    3. Deformable Cylinder lift-off.
    4. Rigid Diamond lift-off.
    ***************/

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args,0,fluids_parameters.COMPRESSIBLE),tests(stream_type_input,data_directory,solid_body_collection),boundary(0),spring_factor(1),rigid_body_collection(solid_body_collection.rigid_body_collection),
        eno_scheme(1),eno_order(2),rk_order(3),cfl_number((T).6),timesplit(false),exact(false),use_slip(false)
    {
        parse_args.Add("-eno_scheme",&eno_scheme,"eno_scheme","eno scheme");
        parse_args.Add("-eno_order",&eno_order,"eno_order","eno order");
        parse_args.Add("-rk_order",&rk_order,"rk_order","runge kutta order");
        parse_args.Add("-cfl",&cfl_number,"CFL","cfl number");
        parse_args.Add("-timesplit",&timesplit,"split time stepping into an explicit advection part, and an implicit non-advection part");
        parse_args.Add("-exact",&exact,"output a fully-explicit sim to (output_dir)_exact");
        parse_args.Add("-spring_factor",&spring_factor,"factor","spring factor");
        parse_args.Add("-slip",&use_slip,"use slip/spd for coupling");
        parse_args.Add_Not("-no_preconditioner",&fluids_parameters.use_preconditioner_for_slip_system,"Do not use preconditioner for slip system");
        parse_args.Parse();

        tests.data_directory=data_directory;
        timesplit=timesplit && !exact;
        if(resolution==1) resolution=50; // Stupid hack for a bad default parameter.

        //grid
        int cells=resolution;
        if(test_number==1 || test_number==3 || test_number==4){
            fluids_parameters.grid->Initialize(TV_INT(5,1)*cells+1,RANGE<VECTOR<T,2> >(VECTOR<T,2>((T)0,(T)0), VECTOR<T,2>((T)1,(T).2)));
            *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid();
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;
    
            if(!user_last_frame) last_frame=350;
if(!this->user_frame_rate) frame_rate=(T)1000.;}
        else{
            fluids_parameters.grid->Initialize(TV_INT(4,3)*cells+1,RANGE<VECTOR<T,2> >(VECTOR<T,2>((T)0,(T)0), VECTOR<T,2>((T)1,(T).75)));
            *fluids_parameters.grid=fluids_parameters.grid->Get_MAC_Grid();
            fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=true;
            fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true;

            if(!user_last_frame) last_frame=10;
            if(!this->user_frame_rate) frame_rate=(T)10000.;}

        fluids_parameters.cfl=cfl_number;
        //custom stuff . . . 
        fluids_parameters.compressible_eos=new EOS_GAMMA<T>;

        if(eno_scheme==1) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,false,false);
        else if(eno_scheme==2) fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,true,false);
        else fluids_parameters.compressible_conservation_method = new CONSERVATION_ENO_LLF<TV,TV::m+2>(true,true,true);
        fluids_parameters.compressible_spatial_order=eno_order;
        //fluids_parameters.compressible_conservation_method->Save_Fluxes();
        //fluids_parameters.compressible_conservation_method->Scale_Outgoing_Fluxes_To_Clamp_Variable(true,0,1e-6);
        fluids_parameters.compressible_rungekutta_order=rk_order;
        fluids_parameters.compressible_timesplit=timesplit;
        fluids_parameters.solve_neumann_regions=false;

        fluids_parameters.use_slip=use_slip;
        fluids_parameters.use_preconditioner_for_slip_system=true;

        if(test_number==1 || test_number==2 || test_number==3 || test_number==4){
            fluids_parameters.solid_affects_fluid=true;
            fluids_parameters.fluid_affects_solid=true;
            solids_fluids_parameters.use_leakproof_solve=false;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            solids_parameters.rigid_body_collision_parameters.use_push_out=true;
            solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T)1;
            solids_parameters.use_trapezoidal_rule_for_velocities=false;
            solids_parameters.verbose_dt=true;
            solid_body_collection.print_residuals=false;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-16;
            solids_parameters.implicit_solve_parameters.cg_projection_iterations=0; // TODO: check this
            solids_parameters.implicit_solve_parameters.cg_iterations=400;
        }

        if(test_number==1 || test_number==3 || test_number==4){
            state_left=TV_DIMENSION((T)5.4,(T)2.2222222222222,(T)0,(T)10.3333333333333);
            state_right=TV_DIMENSION((T)1.4,(T)0,(T)0,(T)1);}
        else if(test_number==2){
            state_left=TV_DIMENSION((T)4.333333333333333,(T)3.2817*(T)sqrt(1e5),(T)0,(T)1.5e6);
            state_right=TV_DIMENSION((T)1,(T)0,(T)0,(T)1.5e5);}


        if(!this->user_output_directory){
            if(timesplit) viewer_dir.output_directory=LOG::sprintf("Standard_Tests/Test_%d__Resolution_%d_%d_semiimplicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
            else viewer_dir.output_directory=LOG::sprintf("Standard_Tests/Test_%d__Resolution_%d_%d_explicit",test_number,(fluids_parameters.grid->counts.x),(fluids_parameters.grid->counts.y));
            if(eno_scheme==2) viewer_dir.output_directory+="_density_weighted";
            else if(eno_scheme==3) viewer_dir.output_directory+="_velocity_weighted";
            if(spring_factor!=(T)1) viewer_dir.output_directory=LOG::sprintf("%s_spring_factor_%lf",viewer_dir.output_directory.c_str(),spring_factor);}

        std::string gnuplot_file=viewer_dir.output_directory+"/common/gnuplot_data.dat";
        gnuplot_file_stream.open(gnuplot_file.c_str());
        LOG::cout<<"writing to file "<<gnuplot_file<<std::endl;
    }

    virtual ~STANDARD_TESTS()
    {
        gnuplot_file_stream.close();
    }

//#####################################################################
// Function Intialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    //set custom boundary
    VECTOR<VECTOR<bool,2>,TV::m> valid_wall;
    for(int axis=0;axis<TV::m;axis++) for(int axis_side=0;axis_side<2;axis_side++)
        valid_wall[axis][axis_side]=(fluids_parameters.mpi_grid?!fluids_parameters.mpi_grid->Neighbor(axis,axis_side):true) && !fluids_parameters.domain_walls[axis][axis_side];

    if(test_number==1 || test_number==3 || test_number==4){
        TV far_field_velocity_left=TV(state_left(1),state_left(2)),far_field_velocity_right=TV(state_right(1),state_right(2));
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,
            T_FACE_VECTOR(state_left(0),state_right(0),state_left(0),state_right(0)),T_FACE_VECTOR(state_left(3),state_right(3),state_left(3),state_right(3)),
            TV_FACE_VECTOR(far_field_velocity_left,far_field_velocity_right,far_field_velocity_left,far_field_velocity_right),(T)1,valid_wall);}
    else if(test_number==2){
        TV far_field_velocity_left=TV(state_left(1),state_left(2)),far_field_velocity_right=TV(state_right(1),state_right(2));
        fluids_parameters.compressible_boundary=new BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>(fluids_parameters.euler,
            T_FACE_VECTOR(state_left(0),state_right(0),state_left(0),state_right(0)),T_FACE_VECTOR(state_left(3),state_right(3),state_left(3),state_right(3)),
            TV_FACE_VECTOR(far_field_velocity_left,far_field_velocity_right,far_field_velocity_left,far_field_velocity_right),(T).5,valid_wall);}
}
//#####################################################################
// Function Intialize_Euler_State
//#####################################################################
void Initialize_Euler_State() override
{
    GRID<TV>& grid=fluids_parameters.euler->grid;
    ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& U=fluids_parameters.euler->U;
    fluids_parameters.euler->e_min=1e-6;

    //initialize grid variables
    for(CELL_ITERATOR<TV> iterator(fluids_parameters.euler->grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        T rho,u_vel,v_vel,p;
        if(test_number==1 || test_number==3 || test_number==4){
            if(grid.X(cell_index).x<=(T).08){rho=state_left(0);u_vel=state_left(1);v_vel=state_left(2);p=state_left(3);}
            else{rho=state_right(0);u_vel=state_right(1);v_vel=state_right(2);p=state_right(3);}}
        else if(test_number==2){
            if(grid.X(cell_index).x<=(T).05){rho=state_left(0);u_vel=state_left(1);v_vel=state_left(2);p=state_left(3);}
            else{rho=state_right(0);u_vel=state_right(1);v_vel=state_right(2);p=state_right(3);}}
        else{rho=(T)1.4;u_vel=(T)0;v_vel=(T)0;p=(T)1;}

        U(cell_index)(0)=rho;U(cell_index)(1)=rho*u_vel;U(cell_index)(2)=rho*v_vel;
        U(cell_index)(3)=rho*(fluids_parameters.euler->eos->e_From_p_And_rho(p,rho)+(sqr(u_vel)+sqr(v_vel))*((T).5));}
}
//#####################################################################
// Function Intialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{   
    if(test_number==1){
        int sphere=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies_2D/circle",(T).05,true,true,false);
        rigid_body_collection.Rigid_Body(sphere).Set_Coefficient_Of_Restitution((T)1);rigid_body_collection.Rigid_Body(sphere).coefficient_of_friction=(T)1;
        rigid_body_collection.Rigid_Body(sphere).Is_Kinematic()=false;rigid_body_collection.rigid_body_particles.frame(sphere).t=TV((T).15,(T).05);
        rigid_body_collection.Rigid_Body(sphere).Set_Mass(pi*.05*.05*10.77);}
    else if(test_number==2){
        TRIANGULATED_AREA<T>& top_boundary=tests.Create_Mattress(GRID<TV>(TV_INT(101,26),RANGE<TV>(TV(0,(T).5),TV((T)1,(T).75))),true,0,(T)3.1538,true);
        SEGMENTED_CURVE_2D<T>& top_boundary_curve=top_boundary.Get_Boundary_Object();
        solid_body_collection.deformable_body_collection.Add_Structure(&top_boundary_curve); // TODO(jontg): Necessary?
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& top_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(top_boundary_curve);
        top_collisions.object.Initialize_Hierarchy();Add_To_Fluid_Simulation(top_collisions);

        TRIANGULATED_AREA<T>& bottom_boundary=tests.Create_Mattress(GRID<TV>(TV_INT(101,26),RANGE<TV>(TV(0,0),TV((T)1,(T).25))),true,0,(T)3.1538);
        SEGMENTED_CURVE_2D<T>& bottom_boundary_curve=bottom_boundary.Get_Boundary_Object();
        solid_body_collection.deformable_body_collection.Add_Structure(&bottom_boundary_curve); // TODO(jontg): Necessary?
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& bottom_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(bottom_boundary_curve);
        bottom_collisions.object.Initialize_Hierarchy();Add_To_Fluid_Simulation(bottom_collisions);

        DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        for(int i=0;i<particles.Size();i++)
            if(particles.X(i).x==(T)0 || particles.X(i).x==(T)1 || particles.X(i).y==(T)0 || (abs(particles.X(i).y-(T).75) < 1e-10)) bound_particles.Append(i);

        original_position=solid_body_collection.deformable_body_collection.particles.X(50);

        solid_body_collection.deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
        solid_body_collection.Add_Force(Create_Edge_Springs(top_boundary,(T)1.5e5*spring_factor,(T)0));solid_body_collection.Add_Force(Create_Edge_Springs(bottom_boundary,(T)1.5e5*spring_factor,(T)0));
        solid_body_collection.Add_Force(Create_Altitude_Springs(top_boundary,(T)1.5e5*spring_factor));solid_body_collection.Add_Force(Create_Altitude_Springs(bottom_boundary,(T)1.5e5*spring_factor));}
    else if(test_number==3){
        TRIANGULATED_AREA<T>& deformable_circle=tests.Create_Triangulated_Object(data_directory+"/Triangulated_Areas/circle-216.tri2d",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T).15,(T).05))),true,false,(T).05);

        // correct number nodes
        for(int i=0;i<solid_body_collection.deformable_body_collection.structures.m;i++){
            solid_body_collection.deformable_body_collection.structures(i)->Update_Number_Nodes();}

        // correct mass
        solid_body_collection.deformable_body_collection.particles.Compute_Auxiliary_Attributes(
            solid_body_collection.deformable_body_collection.soft_bindings);
        
        solid_body_collection.Add_Force(Create_Edge_Springs(deformable_circle,(T).3,(T)1.5));
        solid_body_collection.Add_Force(Create_Altitude_Springs(deformable_circle,(T).3));

        boundary=&deformable_circle.Get_Boundary_Object();

        solid_body_collection.deformable_body_collection.Add_Structure(boundary);
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions=*new DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>(*boundary);
        deformable_collisions.object.Initialize_Hierarchy();
        Add_To_Fluid_Simulation(deformable_collisions);}
    else if(test_number==4){ 
       RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(solid_body_collection.rigid_body_collection,true);
        GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;
        SEGMENTED_CURVE_2D<T>* curve=SEGMENTED_CURVE_2D<T>::Create(particles);

        particles.Add_Elements(8);
        particles.X(0)=TV((T)0     ,(T)-.05);particles.X(1)=TV((T).006 ,(T)-.025);
        particles.X(2)=TV((T).0125 ,(T)0);   particles.X(3)=TV((T).006 ,(T).025);  
        particles.X(4)=TV((T)0     ,(T).05); particles.X(5)=TV((T)-.006,(T).025);  
        particles.X(6)=TV((T)-.0125,(T)0);   particles.X(7)=TV((T)-.006,(T)-.025);

        curve->mesh.number_nodes=particles.Size();curve->mesh.elements.Preallocate(particles.Size());
        for(int i=0;i<particles.Size()-1;i++) curve->mesh.elements.Append(VECTOR<int,2>(i,i+1));
        curve->mesh.elements.Append(VECTOR<int,2>(particles.Size()-1,0));
        curve->Update_Segment_List();curve->Update_Bounding_Box();
        rigid_body.Add_Structure(*curve);

        GRID<TV> *grid=new GRID<TV>(TV_INT(80,80),curve->bounding_box->Thickened((T).025));
        ARRAY<T,VECTOR<int,2> > *phi=new ARRAY<T,VECTOR<int,2> >(grid->Domain_Indices());
        for(CELL_ITERATOR<TV> iterator(*grid);iterator.Valid();iterator.Next())
            (*phi)(iterator.Cell_Index())=curve->Calculate_Signed_Distance(iterator.Location());
        LEVELSET_IMPLICIT_OBJECT<TV> *implicit_object=new LEVELSET_IMPLICIT_OBJECT<TV>(*grid,*phi);
        rigid_body.Add_Structure(*implicit_object);

        solid_body_collection.rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);

        int sphere=rigid_body.particle_index;
        rigid_body.Frame().r=ROTATION<TV>::From_Angle((T)pi/4);
        rigid_body_collection.Rigid_Body(sphere).Set_Coefficient_Of_Restitution((T)1);rigid_body_collection.Rigid_Body(sphere).coefficient_of_friction=(T)1;
        rigid_body_collection.Rigid_Body(sphere).Is_Kinematic()=false;rigid_body_collection.rigid_body_particles.frame(sphere).t=TV((T).15,(T).04);
        
        MASS_PROPERTIES<TV> mass_properties(*curve,true);
        mass_properties.Set_Density((T).1077);
        rigid_body_collection.Rigid_Body(sphere).Mass()=mass_properties.Mass();
        rigid_body_collection.Rigid_Body(sphere).Inertia_Tensor()=DIAGONAL_MATRIX<T,1>(mass_properties.Inertia_Tensor().x00);}

    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);
    if(fluids_parameters.use_slip){
        for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->compute_half_forces=true;
        for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++) solid_body_collection.deformable_body_collection.deformables_forces(k)->compute_half_forces=true;
        for(int i=0;i<solid_body_collection.rigid_body_collection.rigids_forces.m;i++) solid_body_collection.rigid_body_collection.rigids_forces(i)->compute_half_forces=true;}

    THIN_SHELLS_FLUID_COUPLING_UTILITIES<T>::Add_Rigid_Body_Walls(*this);

    // Output, in gnuplot-parsable form, the surface of the mesh
    if(test_number==3){
        std::string gnuplot_file=LOG::sprintf("%s/common/gnuplot_data_0.dat",viewer_dir.output_directory.c_str());
        std::ofstream gnuplot_surface_stream;
        gnuplot_surface_stream.open(gnuplot_file.c_str());
        boundary->mesh.Initialize_Ordered_Loop_Nodes();assert(boundary->mesh.ordered_loop_nodes->m==1);
        ARRAY<int>& segmented_curve=(*boundary->mesh.ordered_loop_nodes)(0);
        for(int i=0;i<segmented_curve.m;i++){
            LOG::cout<<boundary->particles.X(segmented_curve(i))(0)<<"\t"<<boundary->particles.X(segmented_curve(i))(1)<<std::endl;
            gnuplot_surface_stream<<boundary->particles.X(segmented_curve(i))(0)<<"\t"<<boundary->particles.X(segmented_curve(i))(1)<<std::endl;}
        gnuplot_surface_stream.flush();
        gnuplot_surface_stream.close();}
}
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override 
{
    V.Subset(bound_particles).Fill(TV());
}
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override 
{
    V.Subset(bound_particles).Fill(TV());
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
void Postprocess_Substep(const T dt,const T time) override 
{
    T rotation=0;
    TV position(0,0),velocity(0,0);
    if(test_number==1 || test_number==4){
        RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particles;
        int rigid_body_index=0;
        position=rigid_body_particles.frame(rigid_body_index).t;
        velocity=rigid_body_particles.twist(rigid_body_index).linear;
        rotation=rigid_body_particles.twist(rigid_body_index).angular.x;}
    else if(test_number==3){
        rotation=0;
        position=solid_body_collection.deformable_body_collection.particles.Center_Of_Mass();}
    else if(test_number==2){
        position=solid_body_collection.deformable_body_collection.particles.X(50);
        velocity=original_position-solid_body_collection.deformable_body_collection.particles.X(50);
        rotation=velocity.Magnitude();
    }

    LOG::cout<<"writing to file :::::"<<""<<time<<" "<<position<<" "<<velocity<<" "<<rotation<<std::endl;
    gnuplot_file_stream<<""<<time<<" "<<position.x<<" "<<position.y<<" "<<velocity.x<<" "<<velocity.y<<" "<<rotation<<std::endl;
    gnuplot_file_stream.flush();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files() const override
{
    BASE::Write_Output_Files();
    if(test_number==3)
    {
        std::string gnuplot_file=viewer_dir.current_directory+"gnuplot_data.dat";
        std::ofstream gnuplot_surface_stream(gnuplot_file.c_str());
        
        boundary->mesh.Initialize_Ordered_Loop_Nodes();assert(boundary->mesh.ordered_loop_nodes->m==1);
        ARRAY<int>& segmented_curve=(*boundary->mesh.ordered_loop_nodes)(0);
        for(int i=0;i<segmented_curve.m;i++){
            LOG::cout<<boundary->particles.X(segmented_curve(i))(0)<<"\t"<<boundary->particles.X(segmented_curve(i))(1)<<std::endl;
            gnuplot_surface_stream<<boundary->particles.X(segmented_curve(i))(0)<<"\t"<<boundary->particles.X(segmented_curve(i))(1)<<std::endl;}
    }
}
//#####################################################################
};
}
#endif
