//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/TYPED_STREAM.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include "DELAUNAY_TRIANGULATION_2D.h"
#include "MPM_SIMULATION.h"
#include "MPM_PROJECTION.h"
#include "SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL.h"
#include "SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON.h"
#include "TIMING.h"
#include "VORONOI_2D.h"
#include <omp.h>
using namespace PhysBAM;

template<class TV>
void Run_Simulation(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,2> TV_INT;
    typedef float RW;
    MPM_SIMULATION<TV> sim;
    VORONOI_2D<T> voronoi;
    int test_number=-1;
    std::string output_directory="";
    bool use_output_directory=false;
    sim.dt=(T)1e-3;
    int frame_jump=20;
    int grid_res=16,particle_res=32,particle_count=100;
    T particle_exclude_radius=0;
    T density_scale=(T)1;
    T ym=1;
    T pr=(T).3;
    bool use_voronoi=false;
    bool use_voronoi_boundary=false;
    bool use_turk=false;
    bool use_bridson=false;
    bool use_delaunay=false;
    bool use_projection=false;
    T delaunay_maximum_edge_length=(T)99999;
    T delaunay_minimum_angle=(T)0;
    sim.xi=(T)0;
    sim.PROFILING=false;
    parse_args.Add("-test",&test_number,"test","test number");
    parse_args.Add("-o",&output_directory,&use_output_directory,"o","output directory");
    parse_args.Add("-dt",&sim.dt,"dt","dt");
    parse_args.Add("-fj",&frame_jump,"fj","frame jump");
    parse_args.Add("-dump_matrix",&sim.dump_matrix,"dump linear system");
    parse_args.Add("-test_system",&sim.test_system,"test linear system");
    parse_args.Add("-voronoi",&use_voronoi,"use voronoi mesh reconstruction on the fly");
    parse_args.Add("-voronoib",&use_voronoi_boundary,"use voronoi boundary mesh reconstruction on the fly");
    parse_args.Add("-turk",&use_turk,"use turk mesh reconstruction on the fly");
    parse_args.Add("-bridson",&use_bridson,"use bridson mesh reconstruction on the fly");
    parse_args.Add("-delaunay",&use_delaunay,"use delaunay to generate voronoi");
    parse_args.Add("-delaunay_maxl",&delaunay_maximum_edge_length,"value","triangles with edge longer than this will get deleted");
    parse_args.Add("-delaunay_mina",&delaunay_minimum_angle,"value","triangles with angle smaller than this will get deleted");
    parse_args.Add("-stiffness",&ym,"value","scale stiffness");
    parse_args.Add("-poisson_ratio",&pr,"value","poisson's ratio");
    parse_args.Add("-hardening",&sim.xi,"value","harderning coefficient for normal plasticity");
    parse_args.Add("-flip",&sim.FLIP_alpha,"value","flip fraction");
    parse_args.Add("-gres",&grid_res,"value","grid resolution");
    parse_args.Add("-pres",&particle_res,"value","particle resolution");
    parse_args.Add("-pn",&particle_count,"value","particle number");
    parse_args.Add("-exclude",&particle_exclude_radius,"value","particle exclude radius when using pn");
    parse_args.Add("-rho",&density_scale,"value","scale object density");
    parse_args.Add("-profile",&sim.PROFILING,"print out timing statements");
    parse_args.Add("-projection",&use_projection,"use poisson projection to enforce incompressibility");
    parse_args.Parse(true);
    parse_args.Parse();

    if(!use_output_directory) output_directory=STRING_UTILITIES::string_sprintf("MPM_%dD_test_%d",TV::m,test_number);

    // geometry setting
    switch(test_number){
        case 1:{
            sim.grid.Initialize(TV_INT(1.6*grid_res+1,6*grid_res+1),RANGE<TV>(TV(-0.8,-0.6),TV(0.8,5.4)));

            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-0.3,0),TV(-0.1,5)),particle_exclude_radius);
            int first_count=sim.particles.number;
            sim.particles.Set_Material_Properties(0,first_count,
                (T)115.2*density_scale/first_count, // mass per particle
                (3e3)*ym/((T)2*((T)1+pr)), // mu
                (3e3)*ym*pr/(((T)1+pr)*((T)1-2*pr))); // lambda
            sim.particles.Set_Plasticity(0,first_count,
                false,-1,1, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(0,first_count,
                true,700, // visco_nu
                900, // visco_tau
                0); // visco_kappa

            sim.particles.Add_Randomly_Sampled_Object(SPHERE<TV>(TV(0.2,0.4),0.15),particle_exclude_radius);
            int second_count=sim.particles.number-first_count;
            sim.particles.Set_Material_Properties(first_count,second_count,
                (T)115.2*density_scale/second_count, // mass per particle
                (4e3)*ym/((T)2*((T)1+pr)), // mu
                (4e3)*ym*pr/(((T)1+pr)*((T)1-2*pr))); // lambda
            sim.particles.Set_Plasticity(first_count,second_count,
                false,-1,1, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(first_count,second_count,
                false,700, // visco_nu
                900, // visco_tau
                0); // visco_kappa

            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(0.1,0.6),TV(0.3,2)),particle_exclude_radius);
            int third_count=sim.particles.number-first_count-second_count;
            sim.particles.Set_Material_Properties(first_count+second_count,third_count,
                (T)115.2*density_scale/third_count, // mass per particle
                (0)*ym/((T)2*((T)1+pr)), // mu
                (0)*ym*pr/(((T)1+pr)*((T)1-2*pr))); // lambda
            sim.particles.Set_Plasticity(first_count+second_count,third_count,
                false,-1,1, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(first_count+second_count,third_count,
                false,700, // visco_nu
                900, // visco_tau
                0); // visco_kappa

            sim.use_gravity=true;
            sim.particles.Set_Initial_State(0,sim.particles.number,
                MATRIX<T,TV::m>::Identity_Matrix(), // Fe
                MATRIX<T,TV::m>::Identity_Matrix(), // Fp
                TV()); // initial velocity
            break;}

        case 2:{ // intially stretched
            sim.grid.Initialize(TV_INT(1.2*grid_res+1,1.2*grid_res+1),RANGE<TV>(TV(-0.6,-0.6),TV(0.6,0.6)));

            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-0.2,-0.2),TV(0.2,0.2)),particle_exclude_radius);
            int first_count=sim.particles.number;
            sim.particles.Set_Material_Properties(0,first_count,
                (T)115.2*density_scale/first_count, // mass per particle
                (3e3)*ym/((T)2*((T)1+pr)), // mu
                (3e3)*ym*pr/(((T)1+pr)*((T)1-2*pr))); // lambda
            sim.particles.Set_Plasticity(0,first_count,
                false,-1,1, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(0,first_count,
                false,700, // visco_nu
                900, // visco_tau
                0); // visco_kappa
            sim.use_gravity=false;
            sim.particles.Set_Initial_State(0,sim.particles.number,
                MATRIX<T,TV::m>(2,0,0,0.5), // Fe
                MATRIX<T,TV::m>::Identity_Matrix(), // Fp
                TV()); // initial velocity
            for(int p=0;p<sim.particles.number;p++)
                sim.particles.X(p)=MATRIX<T,TV::m>(2,0,0,0.5)*sim.particles.Xm(p);
            break;}

        case 3:{ // free fall
            sim.grid.Initialize(TV_INT(1.2*grid_res+1,1.2*grid_res+1),RANGE<TV>(TV(-0.6,-0.6),TV(0.6,0.6)));

            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-0.2,-0.2),TV(0.2,0.2)),particle_exclude_radius);
            int first_count=sim.particles.number;
            sim.particles.Set_Material_Properties(0,first_count,
                (T)115.2*density_scale/first_count, // mass per particle
                (3e3)*ym/((T)2*((T)1+pr)), // mu
                (3e3)*ym*pr/(((T)1+pr)*((T)1-2*pr))); // lambda
            sim.particles.Set_Plasticity(0,first_count,
                false,-1,1, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(0,first_count,
                false,700, // visco_nu
                900, // visco_tau
                0); // visco_kappa
            sim.use_gravity=true;
            sim.particles.Set_Initial_State(0,sim.particles.number,
                MATRIX<T,TV::m>::Identity_Matrix(), // Fe
                MATRIX<T,TV::m>::Identity_Matrix(), // Fp
                TV()); // initial velocity
            break;}

        default: PHYSBAM_FATAL_ERROR("Missing test");};

    sim.Initialize();
    if(abs(sim.grid.dX(0)-sim.grid.dX(1))>(T)1e-10){
        LOG::cout<<"grid not uniform! dx: "<<sim.grid.dX<<std::endl;
        exit(0);}

    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,output_directory);

    // Delaunay Triangulation
    TRIANGULATED_AREA<T> ta;
    if(use_delaunay){
        DELAUNAY_TRIANGULATION_2D<T>::Triangulate(sim.particles.X,ta,delaunay_maximum_edge_length,delaunay_minimum_angle);
        for(int s=0;s<ta.mesh.elements.m;s++){
            Add_Debug_Object(VECTOR<TV,TV::m>(sim.particles.X(ta.mesh.elements(s)(0)),sim.particles.X(ta.mesh.elements(s)(1))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));
            Add_Debug_Object(VECTOR<TV,TV::m>(sim.particles.X(ta.mesh.elements(s)(1)),sim.particles.X(ta.mesh.elements(s)(2))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));
            Add_Debug_Object(VECTOR<TV,TV::m>(sim.particles.X(ta.mesh.elements(s)(2)),sim.particles.X(ta.mesh.elements(s)(0))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("delaunay triangulation");
        SEGMENTED_CURVE_2D<T>& boundary_curve=ta.Get_Boundary_Object();
        for(int s=0;s<boundary_curve.mesh.elements.m;s++) Add_Debug_Object(VECTOR<TV,TV::m>(ta.particles.X.Subset(boundary_curve.mesh.elements(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("delaunay triangulation boundary");}
    
    // voronoi reconstruction
    if(use_voronoi){
        voronoi.Initialize_With_And_As_A_Triangulated_Area_And_Relocate_Particles_To_Tri_Centers(ta,sim.particles);
        voronoi.Build_Segments();
        for(int i=0;i<voronoi.X.m;i++) {
            if(voronoi.type(i)==1) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,0,0));
            else if(voronoi.type(i)==10) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,1,0));
            else if(voronoi.type(i)==100) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(0,0,1));}
        for(int s=0;s<voronoi.segments.m;s++){
            Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("voronoi cells");}
    if(use_voronoi_boundary){
        voronoi.Initialize_With_And_As_A_Triangulated_Area_And_Relocate_Particles_To_Tri_Centers(ta,sim.particles);
        voronoi.Build_Boundary_Segments();
        for(int i=0;i<voronoi.X.m;i++) {
            if(voronoi.type(i)==1) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,0,0));
            else if(voronoi.type(i)==10) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,1,0));
            else if(voronoi.type(i)==100) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(0,0,1));}
        for(int s=0;s<voronoi.boundary_segments.m;s++){
            Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.boundary_segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("voronoi cell boundary");}

    // projection init
    MPM_PROJECTION<TV> projection(sim);

    // Greg Turk
    if(use_turk){
        SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV> recons;
        ARRAY<TV> xbar;
        ARRAY<MATRIX<T,TV::m> > G;
        ARRAY<T> density;
        GRID<TV> recons_grid(TV_INT(0.6*300+1,2*300+1),RANGE<TV>(TV(-0.3,-1),TV(0.3,1)));
        ARRAY<T,TV_INT> phi;
        recons.Compute_Kernal_Centers_And_Transformation_And_Density(sim.particles.X,sim.particles.mass,0.02,0.04,0.9,10,4,1400,0.5,xbar,G,density);
        recons.Build_Scalar_Field(xbar,sim.particles.mass,density,G,recons_grid,phi);
        T k;std::cin>>k;
        for(int i=0;i<phi.array.m;i++) phi.array(i)-=k;
    }

    // Zhu and Bridson
    SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<TV> recons_bridson;
    GRID<TV> recons_bridson_grid(TV_INT()+(grid_res*8+1),RANGE<TV>(TV(-1,-1),TV(1,1)));
    ARRAY<T,TV_INT> phi_bridson;
    if(use_bridson){
        recons_bridson.Initialize(particle_exclude_radius*2);
        recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,1);
        phi_bridson.array+=1;}

    // draw wall
    for(T x=sim.grid.domain.min_corner(0)+3.5*sim.grid.dX.Min();x<sim.grid.domain.max_corner(0)-3.5*sim.grid.dX.Min();x+=sim.grid.dX.Min()*0.2){
        Add_Debug_Particle(TV(x,sim.grid.domain.min_corner(1)+3.5*sim.grid.dX.Min()),VECTOR<T,3>(0,0,1));
        Add_Debug_Particle(TV(x,sim.grid.domain.max_corner(1)-3.5*sim.grid.dX.Min()),VECTOR<T,3>(0,0,1));}
    for(T y=sim.grid.domain.min_corner(1)+3.5*sim.grid.dX.Min();y<sim.grid.domain.max_corner(1)-3.5*sim.grid.dX.Min();y+=sim.grid.dX.Min()*0.2){
        Add_Debug_Particle(TV(sim.grid.domain.min_corner(0)+3.5*sim.grid.dX.Min(),y),VECTOR<T,3>(0,0,1));
        Add_Debug_Particle(TV(sim.grid.domain.max_corner(0)-3.5*sim.grid.dX.Min(),y),VECTOR<T,3>(0,0,1));}

    // draw collision objects
    for(int b=0;b<sim.rigid_ball.m;b++){
        for(int k=0;k<50;k++){
            T theta=k*2.0*3.14/50.0;
            TV disp;disp(0)=sim.rigid_ball(b).radius*cos(theta);disp(1)=sim.rigid_ball(b).radius*sin(theta);
            Add_Debug_Particle(sim.rigid_ball(b).center+disp,VECTOR<T,3>(0,0,1));}}

    // draw MPM particles
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,1,1));

    Flush_Frame<TV>("mpm");

    for(int f=1;f<2900977;f++){
        TIMING_START;
        LOG::cout<<"MPM TIMESTEP "<<f<<std::endl;

        sim.Build_Weights_And_Grad_Weights();
        sim.Build_Helper_Structures_For_Constitutive_Model();
        LOG::cout<<"Momentum - particles:"<<sim.Get_Total_Momentum_On_Particles()<<std::endl;
        sim.Rasterize_Particle_Data_To_The_Grid();
        LOG::cout<<"Momentum - grid (before linear solve):"<<sim.Get_Total_Momentum_On_Nodes()<<std::endl;
        if(sim.frame==0) sim.Compute_Particle_Volumes_And_Densities();
        sim.Compute_Grid_Forces();
        if(sim.use_gravity) sim.Apply_Gravity_To_Grid_Forces();
        sim.Update_Velocities_On_Grid();
        sim.Grid_Based_Body_Collisions();
        sim.Solve_The_Linear_System(); // so far sim.node_V is achieved via MPM
        LOG::cout<<"Momentum - grid (after linear solve):"<<sim.Get_Total_Momentum_On_Nodes()<<std::endl;
        if(use_projection){
            projection.Reinitialize();
            projection.Identify_Dirichlet_Cells();
            projection.Identify_Neumann_Cells();
            projection.Identify_Nodes_Of_Non_Dirichlet_Cells();
            projection.Velocities_Corners_To_Faces_MPM_Style();
            projection.Build_Velocity_Divergence();
            projection.Solve_For_Pressure(sim.dt,1);
            projection.Do_Projection(sim.dt,1);
            projection.Velocities_Faces_To_Corners_MPM_Style(); // this step modifies sim.node_V
            LOG::cout<<"Momentum - grid (after projection):"<<sim.Get_Total_Momentum_On_Nodes()<<std::endl;
        }

        sim.Update_Deformation_Gradient();
        sim.Update_Particle_Velocities();
        sim.Particle_Based_Body_Collisions();
        sim.Update_Particle_Positions();
        sim.Update_Dirichlet_Box_Positions();
        sim.Update_Colliding_Object_Positions();
        sim.frame++;

        TIMING_END("Current MPM time step totally");

        if(f%frame_jump==0){
            // draw MPM particles
            for(int i=0;i<sim.particles.X.m;i++){
                if(sim.particles.Xm(i).x<0)
                    Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,0,1));
                else if(sim.particles.Xm(i).y<0.55)
                    Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,1,0));
                // else
                    
                    // Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,1));
            }

            // projection: visualize MAC grid velocities
            // if(use_projection){
            //     for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),projection.mac_grid.counts));it.Valid();it.Next()){
            //         Add_Debug_Particle(projection.mac_grid.X(it.index),VECTOR<T,3>(0,0,0)); // cell centers: red
            //         Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.node_V(it.index));}
            //     for(FACE_ITERATOR<TV> iterator(projection.mac_grid);iterator.Valid();iterator.Next()){
            //         TV location=iterator.Location();
            //         int axis=iterator.Axis();
            //         if(axis==0){
            //             Add_Debug_Particle((location),VECTOR<T,3>(0,0,0)); // face centers with x component velocity: green
            //             Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(projection.face_velocities(iterator.Full_Index()),0));}
            //         else if(axis==1){
            //             Add_Debug_Particle((location),VECTOR<T,3>(0,0,0)); // face centers with y component velocity: blue
            //             Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(0,projection.face_velocities(iterator.Full_Index())));}}}

            // Zhu and Bridson
            if(use_bridson){
                recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,1);
                Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(1,0,0));
                LEVELSET<TV> ls(recons_bridson_grid,phi_bridson,0);
                ls.Fast_Marching_Method();
                phi_bridson.array+=particle_exclude_radius*3;
                ls.Fast_Marching_Method();
                phi_bridson.array-=particle_exclude_radius*2;
                ls.Fast_Marching_Method();
                Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(0,1,0));}

            // voronoi reconstruction
            if(use_voronoi){
                voronoi.Crack(sim.particles.X,sim.grid.dX.Min()*1.8);
                voronoi.Build_Association();
                voronoi.Build_Segments();
                voronoi.Deform_Mesh_Using_Particle_Deformation(sim.particles.Xm,sim.particles.X,sim.particles.Fe,sim.particles.Fp,true,3);
                for(int s=0;s<voronoi.segments.m;s++){
                    Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}}
            if(use_voronoi_boundary){
                voronoi.Crack(sim.particles.X,sim.grid.dX.Min()*1.3);
                voronoi.Build_Association();
                voronoi.Build_Boundary_Segments();
                voronoi.Deform_Mesh_Using_Particle_Deformation(sim.particles.Xm,sim.particles.X,sim.particles.Fe,sim.particles.Fp,true,1);
                for(int s=0;s<voronoi.boundary_segments.m;s++){
                    Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.boundary_segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}}

            // draw wall
            for(T x=sim.grid.domain.min_corner(0)+4.5*sim.grid.dX.Min();x<sim.grid.domain.max_corner(0)-4.5*sim.grid.dX.Min();x+=sim.grid.dX.Min()*0.2){
                Add_Debug_Particle(TV(x,sim.grid.domain.min_corner(1)+4.5*sim.grid.dX.Min()),VECTOR<T,3>(0,0,1));
                Add_Debug_Particle(TV(x,sim.grid.domain.max_corner(1)-4.5*sim.grid.dX.Min()),VECTOR<T,3>(0,0,1));}
            for(T y=sim.grid.domain.min_corner(1)+4.5*sim.grid.dX.Min();y<sim.grid.domain.max_corner(1)-4.5*sim.grid.dX.Min();y+=sim.grid.dX.Min()*0.2){
                Add_Debug_Particle(TV(sim.grid.domain.min_corner(0)+4.5*sim.grid.dX.Min(),y),VECTOR<T,3>(0,0,1));
                Add_Debug_Particle(TV(sim.grid.domain.max_corner(0)-4.5*sim.grid.dX.Min(),y),VECTOR<T,3>(0,0,1));}
            
            // draw collision objects
            for(int b=0;b<sim.rigid_ball.m;b++){
                for(int k=0;k<50;k++){
                    T theta=k*2.0*3.14/50.0;
                    TV disp;disp(0)=sim.rigid_ball(b).radius*cos(theta);disp(1)=sim.rigid_ball(b).radius*sin(theta);
                    Add_Debug_Particle(sim.rigid_ball(b).center+disp,VECTOR<T,3>(0,0,1));}}

            Flush_Frame<TV>("mpm");}
        LOG::cout<<std::endl;}
}

int main(int argc,char *argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse(true);
    Run_Simulation<VECTOR<double,2> >(parse_args);
    LOG::Finish_Logging();
    return 0;
}


