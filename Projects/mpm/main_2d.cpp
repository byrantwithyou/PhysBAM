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
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include "DELAUNAY_TRIANGULATION_2D.h"
#include "MPM_PROJECTION.h"
#include "MPM_SIMULATION.h"
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
    parse_args.Add("-projection",&use_projection,"use poisson projection to enforce incompressibility");
    parse_args.Parse(true);
    parse_args.Parse();

    if(!use_output_directory) output_directory=STRING_UTILITIES::string_sprintf("MPM_%dD_test_%d",TV::m,test_number);

    // geometry setting
    switch(test_number){
        case 1:{ // all materials together
            sim.grid.Initialize(TV_INT(2*grid_res+1,2*grid_res+1),RANGE<TV>(TV(-1,-1),TV(1,1)));
            
            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-0.2,-0.2),TV(0.2,0.2)),particle_exclude_radius);
            
            int c1=sim.particles.number;
            sim.particles.Set_Material_Properties(0,c1,
                (T)8*density_scale/c1, // mass per particle
                80.0*ym/(2.0*(1.0+pr)), // mu
                80.0*ym*pr/((1.0+pr)*(1.0-2.0*pr)), // lambda
                true,0); // compress, pressure
            
//           sim.particles.Set_Material_Properties(0,c1,
//                                                 (T)8*density_scale/1000, // mass per particle
//                                                 0, // mu
//                                                 0, // lambda
//                                                 false,0); // compress, pressure
            
            sim.particles.Set_Plasticity(0,c1,
                false,-1000,1.2, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(0,c1,
                false,100, // visco_nu
                5000, // visco_tau
                0); // visco_kappa
            sim.particles.Set_Initial_State(0,c1,
                MATRIX<T,TV::m>(1.5,0,0,0.8), // Fe
                MATRIX<T,TV::m>::Identity_Matrix(), // Fp
                TV(0,0)); // initial velocity
            for(int p=0;p<c1;p++) sim.particles.X(p)=MATRIX<T,TV::m>(1.5,0,0,0.8)*sim.particles.Xm(p);
            
            sim.use_gravity=false;

            break;}
     
        default: PHYSBAM_FATAL_ERROR("Missing test");};

    sim.Initialize();
    if(abs(sim.grid.dX(0)-sim.grid.dX(1))>(T)1e-15){
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

        int c1=sim.particles.number;
        sim.particles.Set_Plasticity(0,c1,
            true,-1000,1.2, // plasticity_yield
            false,-1,1); // plasticity_clamp
        sim.particles.Set_Visco_Plasticity(0,c1,
            false,100, // visco_nu
            5000, // visco_tau
            0); // visco_kappa
        sim.particles.Set_Initial_State(0,c1,
            MATRIX<T,TV::m>::Identity_Matrix(), // Fe
            MATRIX<T,TV::m>::Identity_Matrix(), // Fp
            TV()); // initial velocity
        for(int p=0;p<sim.particles.number;p++){
            T fullE=3000*ym;
            T alpha=1;
            T thisE=(1.0-alpha)*fullE/(0.25*0.25)*sqr(sim.particles.Xm(p).y)+alpha*fullE;
            sim.particles.Set_Material_Properties(p,1,
                 (T)200*density_scale/1000, // mass per particle
                thisE/((T)2*((T)1+pr)), // mu
                thisE*pr/(((T)1+pr)*((T)1-2*pr)), // lambda
                false,0); // compress,pressure
        }
        sim.dirichlet_box.Append(RANGE<TV>(TV(-0.3,-0.3),TV(0.3,-0.22)));
        sim.dirichlet_velocity.Append(TV(0,-0.1));
        sim.dirichlet_box.Append(RANGE<TV>(TV(-0.3,0.22),TV(0.3,0.3)));
        sim.dirichlet_velocity.Append(TV(0,0.1));

        sim.Initialize();


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
    GRID<TV> recons_bridson_grid(TV_INT()+(grid_res*8+1),RANGE<TV>(TV(-0.6,-0.6),TV(0.6,0.6)));
    ARRAY<T,TV_INT> phi_bridson;
    if(use_bridson){
        recons_bridson.Initialize(particle_exclude_radius*4);
        recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,0);
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

    for(int f=0;f<=20032323;f++){

        TIMING_START;
        LOG::cout<<"MPM TIMESTEP "<<f<<std::endl;
        sim.Build_Weights_And_Grad_Weights();
        sim.Build_Helper_Structures_For_Constitutive_Model();
        sim.Build_Pressure_And_One_Over_Lambda_J();
        LOG::cout<<"Momentum - particles:"<<sim.Get_Total_Momentum_On_Particles()<<std::endl;

        // LOG::cout<<"particle_velocity" <<sim.particles.V<<std::endl;

        sim.Rasterize_Particle_Data_To_The_Grid();
        
        // draw particles
//        for(int i=0;i<sim.particles.X.m;i++){
//            Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,1));
//            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.particles.V(i));
//        }
//        // draw node velocity
//        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),sim.grid.counts));it.Valid();it.Next()){
//            Add_Debug_Particle(sim.grid.Node(it.index),VECTOR<T,3>(1,1,1));
//            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.node_V(it.index));
//        }
//        std::string title=STRING_UTILITIES::string_sprintf("TIME STEP %d Before MPM Solve",f);
//        Flush_Frame<TV>(title.c_str());

        // LOG::cout<<"Momentum - grid (before linear solve):"<<sim.Get_Total_Momentum_On_Nodes()<<std::endl;
        // LOG::cout<<"node velocity"<<sim.node_V<<std::endl;

        if(sim.frame==0) sim.Compute_Particle_Volumes_And_Densities(0,sim.particles.number);

        sim.Compute_Grid_Forces();
        // sim.node_force.Fill(TV());
        
        if(sim.use_gravity) sim.Apply_Gravity_To_Grid_Forces();
        sim.Update_Velocities_On_Grid();
        if(!use_projection) sim.Grid_Based_Body_Collisions();

        sim.Solve_The_Linear_System(); // so far sim.node_V is achieved via MPM
        // LOG::cout<<"Momentum - grid (after linear solve):"<<sim.Get_Total_Momentum_On_Nodes()<<std::endl;

//        // draw particles
//        for(int i=0;i<sim.particles.X.m;i++){
//            Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,1));
//            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.particles.V(i));
//        }
//        // draw node velocity
//        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),sim.grid.counts));it.Valid();it.Next()){
//            Add_Debug_Particle(sim.grid.Node(it.index),VECTOR<T,3>(1,1,1));
//            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.node_V(it.index));
//        }
//        title=STRING_UTILITIES::string_sprintf("TIME STEP %d After MPM Solve",f);
//        Flush_Frame<TV>(title.c_str());


        if(use_projection){
            projection.Reinitialize();
            projection.Identify_Dirichlet_Cells();
            
            projection.Identify_Neumann_Cells();
            for(int x=3;x<projection.mac_grid.counts.x-3;x++){
                for(int y=2;y<=4;y++){
                    projection.cell_neumann(TV_INT(x,y))=true;
                    projection.cell_dirichlet(TV_INT(x,y))=false;
                    projection.neumann_cell_normal_axis(TV_INT(x,y))=2;}}
            for(int x=3;x<projection.mac_grid.counts.x-3;x++){
                for(int y=projection.mac_grid.counts.y-5;y<=projection.mac_grid.counts.y-3;y++){
                    projection.cell_neumann(TV_INT(x,y))=true;
                    projection.cell_dirichlet(TV_INT(x,y))=false;
                    projection.neumann_cell_normal_axis(TV_INT(x,y))=-2;}}
            for(int x=2;x<=4;x++){
                for(int y=3;y<projection.mac_grid.counts.y-3;y++){
                    projection.cell_neumann(TV_INT(x,y))=true;
                    projection.cell_dirichlet(TV_INT(x,y))=false;
                    projection.neumann_cell_normal_axis(TV_INT(x,y))=1;}}
            for(int x=projection.mac_grid.counts.x-5;x<=projection.mac_grid.counts.x-3;x++){
                for(int y=3;y<projection.mac_grid.counts.y-3;y++){
                    projection.cell_neumann(TV_INT(x,y))=true;
                    projection.cell_dirichlet(TV_INT(x,y))=false;
                    projection.neumann_cell_normal_axis(TV_INT(x,y))=-1;}}
            
            projection.Identify_Nodes_Of_Non_Dirichlet_Cells();
            projection.Velocities_Corners_To_Faces_MPM_Style();


            // LOG::cout<<"face velocity from node"<<projection.face_velocities<<std::endl;

//            // draw particles
//            for(int i=0;i<sim.particles.X.m;i++){
//                Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,1));
//                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.particles.V(i));
//            }
//            // draw node velocity
//            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),sim.grid.counts));it.Valid();it.Next()){
//                Add_Debug_Particle(sim.grid.Node(it.index),VECTOR<T,3>(1,1,1));
//                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.node_V(it.index));
//            }
//            // draw mac velocity
//            for(FACE_ITERATOR<TV> iterator(projection.mac_grid);iterator.Valid();iterator.Next()){
//                TV location=iterator.Location();
//                int axis=iterator.Axis();
//                if(axis==0){
//                    Add_Debug_Particle((location),VECTOR<T,3>(1,1,1)); // face centers with x component velocity: 
//                    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(projection.face_velocities(iterator.Full_Index()),0));}
//                else if(axis==1){
//                    Add_Debug_Particle((location),VECTOR<T,3>(1,1,1)); // face centers with y component velocity: 
//                    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(0,projection.face_velocities(iterator.Full_Index())));}}
//            // visualize Neumann cells and Dirichlet cells
//            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),projection.mac_grid.counts));it.Valid();it.Next()){
//                if(projection.cell_dirichlet(it.index))
//                    Add_Debug_Particle(projection.mac_grid.Center(it.index),VECTOR<T,3>(0.3,0.3,0.3));
//                if(projection.cell_neumann(it.index)){
//                    PHYSBAM_ASSERT(!projection.cell_dirichlet(it.index));
//                    Add_Debug_Particle(projection.mac_grid.Center(it.index),VECTOR<T,3>(0,1,0));}}
//            title=STRING_UTILITIES::string_sprintf("TIME STEP %d After Velocity Corner->Face",f);
//            Flush_Frame<TV>(title.c_str());


            projection.Build_Weights_And_Grad_Weights_For_Cell_Centers();
            projection.Rasterize_Pressure_And_One_Over_Lambda_J();
            projection.Build_Velocity_Divergence();
            LOG::cout<<"Maximum velocity divergence before projection: "<<projection.max_div<<std::endl;
            // LOG::cout<<"div_u"<<projection.div_u<<std::endl;

            projection.Solve_For_Pressure();
            projection.Do_Projection();
            projection.Build_Velocity_Divergence();
            LOG::cout<<"Maximum velocity divergence after projection: "<<projection.max_div<<std::endl;

//            // draw particles
//            for(int i=0;i<sim.particles.X.m;i++){
//                Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,1));
//                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.particles.V(i));
//            }
//            // draw node velocity
//            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),sim.grid.counts));it.Valid();it.Next()){
//                Add_Debug_Particle(sim.grid.Node(it.index),VECTOR<T,3>(1,1,1));
//                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.node_V(it.index));
//            }
//            // draw mac velocity
//            for(FACE_ITERATOR<TV> iterator(projection.mac_grid);iterator.Valid();iterator.Next()){
//                TV location=iterator.Location();
//                int axis=iterator.Axis();
//                if(axis==0){
//                    Add_Debug_Particle((location),VECTOR<T,3>(1,1,1)); // face centers with x component velocity: 
//                    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(projection.face_velocities(iterator.Full_Index()),0));}
//                else if(axis==1){
//                    Add_Debug_Particle((location),VECTOR<T,3>(1,1,1)); // face centers with y component velocity: 
//                    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(0,projection.face_velocities(iterator.Full_Index())));}}
//            // visualize Neumann cells and Dirichlet cells
//            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),projection.mac_grid.counts));it.Valid();it.Next()){
//                if(projection.cell_dirichlet(it.index))
//                    Add_Debug_Particle(projection.mac_grid.Center(it.index),VECTOR<T,3>(0.3,0.3,0.3));
//                if(projection.cell_neumann(it.index)){
//                    PHYSBAM_ASSERT(!projection.cell_dirichlet(it.index));
//                    Add_Debug_Particle(projection.mac_grid.Center(it.index),VECTOR<T,3>(0,1,0));}}
//            title=STRING_UTILITIES::string_sprintf("TIME STEP %d After Projection",f);
//            Flush_Frame<TV>(title.c_str());

            projection.Velocities_Faces_To_Corners_MPM_Style((T)0.95); // this step modifies sim.node_V

//            // draw particles
//            for(int i=0;i<sim.particles.X.m;i++){
//                Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,1));
//                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.particles.V(i));
//            }
//            // draw node velocity
//            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),sim.grid.counts));it.Valid();it.Next()){
//                Add_Debug_Particle(sim.grid.Node(it.index),VECTOR<T,3>(1,1,1));
//                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.node_V(it.index));
//            }
//            // draw mac velocity
//            for(FACE_ITERATOR<TV> iterator(projection.mac_grid);iterator.Valid();iterator.Next()){
//                TV location=iterator.Location();
//                int axis=iterator.Axis();
//                if(axis==0){
//                    Add_Debug_Particle((location),VECTOR<T,3>(1,1,1)); // face centers with x component velocity: 
//                    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(projection.face_velocities(iterator.Full_Index()),0));}
//                else if(axis==1){
//                    Add_Debug_Particle((location),VECTOR<T,3>(1,1,1)); // face centers with y component velocity: 
//                    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,TV(0,projection.face_velocities(iterator.Full_Index())));}}
//            // visualize Neumann cells and Dirichlet cells
//            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),projection.mac_grid.counts));it.Valid();it.Next()){
//                if(projection.cell_dirichlet(it.index))
//                    Add_Debug_Particle(projection.mac_grid.Center(it.index),VECTOR<T,3>(0.3,0.3,0.3));
//                if(projection.cell_neumann(it.index)){
//                    PHYSBAM_ASSERT(!projection.cell_dirichlet(it.index));
//                    Add_Debug_Particle(projection.mac_grid.Center(it.index),VECTOR<T,3>(0,1,0));}}
//            title=STRING_UTILITIES::string_sprintf("TIME STEP %d After Velocity Face To Corner",f);
//            Flush_Frame<TV>(title.c_str());

            projection.Pressure_Back_To_Particles((T)0.95);
            
            // LOG::cout<<"particle pressure:"<<sim.particles.pressure<<std::endl;
            // LOG::cout<<"particle one_over_lambda_J:"<<sim.particles.one_over_lambda_J<<std::endl;

        }

        sim.Update_Deformation_Gradient();
        // LOG::cout<<sim.particles.Fe<<std::endl;

        sim.Update_Particle_Velocities();
        if(!use_projection) sim.Particle_Based_Body_Collisions();
        sim.Update_Particle_Positions();

        sim.Update_Dirichlet_Box_Positions();
        sim.Update_Colliding_Object_Positions();
        sim.frame++;

        TIMING_END("Current MPM time step totally");

        if(f%frame_jump==0){

            // draw MPM particles
            for(int i=0;i<sim.particles.X.m;i++){
                Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,1));}

            // visualize Neumann cells and Dirichlet cells
            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),projection.mac_grid.counts));it.Valid();it.Next()){
        //         if(projection.cell_dirichlet(it.index))
        //             Add_Debug_Particle(projection.mac_grid.Center(it.index),VECTOR<T,3>(0.3,0.3,0.3));
                if(projection.cell_neumann(it.index)){
                    PHYSBAM_ASSERT(!projection.cell_dirichlet(it.index));
                    Add_Debug_Particle(projection.mac_grid.Center(it.index),VECTOR<T,3>(0,1,0));}}

        //     // Zhu and Bridson
        //     if(use_bridson){
        //         recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,0);
        //         Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(1,0,0));
        //         LEVELSET<TV> ls(recons_bridson_grid,phi_bridson,0);
        //         ls.Fast_Marching_Method();
        //         phi_bridson.array+=particle_exclude_radius*3;
        //         ls.Fast_Marching_Method();
        //         phi_bridson.array-=particle_exclude_radius*2;
        //         ls.Fast_Marching_Method();
        //         Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(0,1,0));}

        //     // voronoi reconstruction
        //     if(use_voronoi){
        //         voronoi.Crack(sim.particles.X,sim.grid.dX.Min()*1.2);
        //         voronoi.Build_Association();
        //         voronoi.Build_Segments();
        //         voronoi.Deform_Mesh_Using_Particle_Deformation(sim.particles.Xm,sim.particles.X,sim.particles.Fe,sim.particles.Fp,true,3);
        //         for(int s=0;s<voronoi.segments.m;s++){
        //             Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}}
        //     if(use_voronoi_boundary){
        //         voronoi.Crack(sim.particles.X,sim.grid.dX.Min()*1.3);
        //         voronoi.Build_Association();
        //         voronoi.Build_Boundary_Segments();
        //         voronoi.Deform_Mesh_Using_Particle_Deformation(sim.particles.Xm,sim.particles.X,sim.particles.Fe,sim.particles.Fp,true,1);
        //         for(int s=0;s<voronoi.boundary_segments.m;s++){
        //             Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.boundary_segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}}
            
        //     // draw collision objects
        //     for(int b=0;b<sim.rigid_ball.m;b++){
        //         for(int k=0;k<50;k++){
        //             T theta=k*2.0*3.14/50.0;
        //             TV disp;disp(0)=sim.rigid_ball(b).radius*cos(theta);disp(1)=sim.rigid_ball(b).radius*sin(theta);
        //             Add_Debug_Particle(sim.rigid_ball(b).center+disp,VECTOR<T,3>(0,0,1));}}

        std::string title=STRING_UTILITIES::string_sprintf("TIME STEP %d After Moving Particles",f);
        Flush_Frame<TV>(title.c_str());}

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


