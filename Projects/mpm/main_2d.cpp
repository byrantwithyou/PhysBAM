//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Tools/Read_Write/TYPED_STREAM.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <omp.h>
#include "TIMING.h"
#include "VORONOI_2D.h"
#include "MPM_SIMULATION.h"
#include "SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL.h"
#include "SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON.h"
#include "DELAUNAY_TRIANGULATION_2D.h"
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
    T delaunay_maximum_edge_length=(T)99999;
    T delaunay_minimum_angle=(T)0;
    sim.xi=(T)0;
    sim.visco_kappa=(T)0;
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
    parse_args.Add("-visco_hardening",&sim.visco_kappa,"value","harderning coefficient for visco plasticity");
    parse_args.Add("-flip",&sim.FLIP_alpha,"value","flip fraction");
    parse_args.Add("-gres",&grid_res,"value","grid resolution");
    parse_args.Add("-pres",&particle_res,"value","particle resolution");
    parse_args.Add("-pn",&particle_count,"value","particle number");
    parse_args.Add("-exclude",&particle_exclude_radius,"value","particle exclude radius when using pn");
    parse_args.Add("-rho",&density_scale,"value","scale object density");
    parse_args.Add("-profile",&sim.PROFILING,"print out timing statements");
    parse_args.Parse(true);
    parse_args.Parse();

    if(!use_output_directory) output_directory=STRING_UTILITIES::string_sprintf("MPM_%dD_test_%d",TV::m,test_number);

    T object_mass=1;

    // geometry setting
    switch(test_number){
        case 1: // stretching beam
            sim.grid.Initialize(TV_INT(2*grid_res+1,0.5*grid_res+1),RANGE<TV>(TV(-1,-0.25),TV(1,0.25)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.4*particle_res+1,0.24*particle_res+1),RANGE<TV>(TV(-0.2,-0.12),TV(0.2,0.12)));
            sim.ground_level=-100;
            sim.dirichlet_box.Append(RANGE<TV>(TV(0.16,-0.15),TV(0.25,0.15)));
            sim.dirichlet_velocity.Append(TV(0.2,0));
            sim.dirichlet_box.Append(RANGE<TV>(TV(-0.25,-0.15),TV(-0.16,0.15)));
            sim.dirichlet_velocity.Append(TV(-0.2,0));
            break;
        case 2: // shit fall
            sim.grid.Initialize(TV_INT(1.2*grid_res+1,1.2*grid_res+1),RANGE<TV>(TV(-0.6,-0.6),TV(0.6,0.6)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.2*particle_res+1,0.2*particle_res+1),RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1)));
            sim.particles.Add_X_As_A_Grid(TV_INT(0.2*particle_res+1,0.2*particle_res+1),RANGE<TV>(TV(-0.1,0.2),TV(0.1,0.4)));
            sim.ground_level=-0.55;
            sim.rigid_ball.Append(SPHERE<TV>(TV(0,-0.2),0.03));
            sim.rigid_ball_velocity.Append(TV());
            sim.rigid_ball.Append(SPHERE<TV>(TV(0.09,-0.3),0.03));
            sim.rigid_ball_velocity.Append(TV());
            break;
        case 3: // snow fall
            sim.grid.Initialize(TV_INT(1*grid_res+1,1*grid_res+1),RANGE<TV>(TV(-0.5,-0.7),TV(0.5,0.3)));
            sim.particles.Initialize_X_As_A_Randomly_Sampled_Box(particle_count,RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1)));
            sim.ground_level=-0.6;
            sim.rigid_ball.Append(SPHERE<TV>(TV(0,-0.3),0.05));
            sim.rigid_ball_velocity.Append(TV());
            sim.rigid_ball.Append(SPHERE<TV>(TV(0.13,-0.3),0.05));
            sim.rigid_ball_velocity.Append(TV());
            break;
        case 4: // notch
            sim.grid.Initialize(TV_INT(2*grid_res+1,2*grid_res+1),RANGE<TV>(TV(-1,-1),TV(1,1)));
            // sim.particles.Initialize_X_As_A_Grid(TV_INT(0.2*particle_res+1,0.4*particle_res+1),RANGE<TV>(TV(-0.1,-0.2),TV(0.1,0.2)));
            sim.particles.Initialize_X_As_A_Randomly_Sampled_Box(particle_count,RANGE<TV>(TV(-0.1,-0.2),TV(0.1,0.2)),particle_exclude_radius);
            // sim.particles.Initialize_X_As_A_Randomly_Sampled_Box(particle_count,RANGE<TV>(TV(-0.1,-0.2),TV(0.1,0.2)));
            sim.particles.Reduce_X_As_A_Ball(RANGE<TV>(TV(0.06,-0.04),TV(0.14,0.04)));
            sim.ground_level=-100;
            sim.dirichlet_box.Append(RANGE<TV>(TV(-0.15,-0.23),TV(0.15,-0.17)));
            sim.dirichlet_velocity.Append(TV(0,-0.2));
            sim.dirichlet_box.Append(RANGE<TV>(TV(-0.15,0.17),TV(0.15,0.23)));
            sim.dirichlet_velocity.Append(TV(0,0.2));
            break;
        case 5: // wall
            sim.grid.Initialize(TV_INT(6*grid_res+1,1*grid_res+1),RANGE<TV>(TV(-0.1,-0.2),TV(5.9,0.8)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.1*particle_res+1,0.5*particle_res+1),RANGE<TV>(TV(0,0),TV(0.1,0.5)));
            sim.ground_level=0;
            sim.dirichlet_box.Append(RANGE<TV>(TV(-0.02,-0.02),TV(0.12,0.02)));
            sim.dirichlet_velocity.Append(TV());
            sim.rigid_ball.Append(SPHERE<TV>(TV(-0.1,0.25),0.05));
            sim.rigid_ball_velocity.Append(TV(50,0));
            break;
        default: PHYSBAM_FATAL_ERROR("Missing test");};

    // material setting
    sim.Initialize();
    if(abs(sim.grid.dX(0)-sim.grid.dX(1))>(T)1e-10){
        LOG::cout<<"grid dx: "<<sim.grid.dX<<std::endl;
        exit(0);}
    switch(test_number){
        case 1:
            object_mass=(T)1200*density_scale*RANGE<TV>(TV(-0.2,-0.12),TV(0.2,0.12)).Size();
            ym*=(T)1e3;
            for(int p=0;p<sim.particles.number;p++){
                T this_ym=0.7*ym/0.04*sqr(sim.particles.X(p)(0))+0.3*ym;
                sim.mu(p)=(this_ym/((T)2*((T)1+pr)));
                sim.lambda(p)=(this_ym*pr/(((T)1+pr)*((T)1-2*pr)));}
            sim.use_gravity=false;
            sim.use_plasticity_yield=true;
            sim.yield_min=-100;
            sim.yield_max=1.3;
            break;
        case 2:
            object_mass=(T)1200*density_scale*(RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1))).Size()*2;
            ym*=(T)1e5;
            sim.mu.Fill(ym/((T)2*((T)1+pr)));
            sim.lambda.Fill(ym*pr/(((T)1+pr)*((T)1-2*pr)));
            sim.use_visco_plasticity=true;
            sim.visco_nu=1e3;
            sim.visco_tau.Fill(1000);
            break;
        case 3:
            object_mass=(T)400*density_scale*(RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1))).Size();
            ym*=(T)140000;
            sim.mu.Fill(ym/((T)2*((T)1+pr)));
            sim.lambda.Fill(ym*pr/(((T)1+pr)*((T)1-2*pr)));
            sim.use_plasticity_yield=true;
            sim.yield_min=(T)1-(T)0.025;
            sim.yield_max=(T)1+(T)0.0075;
            break;
        case 4:
            object_mass=(T)1200*density_scale*(RANGE<TV>(TV(-0.1,-0.2),TV(0.1,0.2)).Size()-SPHERE<TV>(TV(0.1,0),0.04).Size());
            ym*=(T)5e3;
            sim.mu.Fill(ym/((T)2*((T)1+pr)));
            sim.lambda.Fill(ym*pr/(((T)1+pr)*((T)1-2*pr)));
            sim.use_gravity=false;
            sim.use_plasticity_yield=true;
            sim.yield_min=-100;
            sim.yield_max=1.5;
            break;
        case 5:
            object_mass=(T)1200*density_scale*(RANGE<TV>(TV(0,0),TV(0.1,0.5)).Size());
            ym*=(T)5e3;
            sim.mu.Fill(ym/((T)2*((T)1+pr)));
            sim.lambda.Fill(ym*pr/(((T)1+pr)*((T)1-2*pr)));
            sim.use_gravity=true;
            sim.use_plasticity_yield=true;
            sim.yield_min=0.8;
            sim.yield_max=1.2;
            break;
        default: PHYSBAM_FATAL_ERROR("Missing test");};
    for(int p=0;p<sim.particles.number;p++){
        sim.particles.mass(p)=object_mass/sim.particles.number;
        sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
        sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
    sim.mu0=sim.mu;
    sim.lambda0=sim.lambda;

    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,output_directory);

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
        Dump_Levelset(recons_grid,phi,VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));}

    // Zhu and Bridson
    SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<TV> recons_bridson;
    GRID<TV> recons_bridson_grid(TV_INT(2*1200+1,2*1200+1),RANGE<TV>(TV(-1,-1),TV(1,1)));
    ARRAY<T,TV_INT> phi_bridson;
    if(use_bridson){
        recons_bridson.Initialize((T)1/(T)particle_res);
        recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,1);
        Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));}

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
        voronoi.Initialize_With_A_Triangulated_Area(ta);
        voronoi.Build_Segments();
        for(int s=0;s<voronoi.segments.m;s++){
            Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("voronoi cells");}
    if(use_voronoi_boundary){
        voronoi.Initialize_With_A_Triangulated_Area(ta);
        voronoi.Build_Boundary_Segments();
        for(int s=0;s<voronoi.boundary_segments.m;s++){
            Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.boundary_segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}
        for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
        Flush_Frame<TV>("voronoi cell boundary");}

    // collision objects
    for(int b=0;b<sim.rigid_ball.m;b++){
        for(int k=0;k<50;k++){
            T theta=k*2.0*3.14/50.0;
            TV disp;disp(0)=sim.rigid_ball(b).radius*cos(theta);disp(1)=sim.rigid_ball(b).radius*sin(theta);
            Add_Debug_Particle(sim.rigid_ball(b).center+disp,VECTOR<T,3>(0,0,1));}}
    for(int i=0;i<50;i++) Add_Debug_Particle(TV(-1+i*0.04,sim.ground_level),VECTOR<T,3>(0,0,1));

    // MPM particles
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));

    Flush_Frame<TV>("mpm");

    for(int f=1;f<2900977;f++){
        TIMING_START;
        LOG::cout<<"TIMESTEP "<<f<<std::endl;
        sim.Advance_One_Time_Step_Backward_Euler();
        TIMING_END("Current time step totally");
        if(f%frame_jump==0){
            // MPM particles
            for(int i=0;i<sim.particles.X.m;i++) if(sim.valid(i)) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));

            // Zhu and Bridson
            if(use_bridson){
                recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,1);
                Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));}

            // voronoi reconstruction
            if(use_voronoi){
                voronoi.Crack(sim.particles.X,sim.grid.dX.Min()*1.3);
                voronoi.Build_Association();
                voronoi.Build_Segments();
                voronoi.Deform_Mesh_Using_Particle_Deformation(sim.particles.Xm,sim.particles.X,sim.particles.Fe,sim.particles.Fp);
                for(int s=0;s<voronoi.segments.m;s++){
                    Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}}
            if(use_voronoi_boundary){
                voronoi.Crack(sim.particles.X,sim.grid.dX.Min()*1.3);
                voronoi.Build_Association();
                voronoi.Build_Boundary_Segments();
                voronoi.Deform_Mesh_Using_Particle_Deformation(sim.particles.Xm,sim.particles.X,sim.particles.Fe,sim.particles.Fp);
                for(int s=0;s<voronoi.boundary_segments.m;s++){
                    Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.boundary_segments(s))),VECTOR<T,3>(1,0.57,0.25),VECTOR<T,3>(0,0,0));}}
            
            // collision objects
            for(int b=0;b<sim.rigid_ball.m;b++){
                for(int k=0;k<50;k++){
                    T theta=k*2.0*3.14/50.0;
                    TV disp;disp(0)=sim.rigid_ball(b).radius*cos(theta);disp(1)=sim.rigid_ball(b).radius*sin(theta);
                    Add_Debug_Particle(sim.rigid_ball(b).center+disp,VECTOR<T,3>(0,0,1));}}
            for(int i=0;i<50;i++) Add_Debug_Particle(TV(-1+i*0.04,sim.ground_level),VECTOR<T,3>(0,0,1));

            Flush_Frame<TV>("mpm");}
        LOG::cout<<std::endl;}
}

int main(int argc,char *argv[])
{
    // TRIANGULATED_AREA<double> ta;
    // ta.Clean_Memory();
    // ta.particles.Delete_All_Elements();
    // ta.particles.Add_Elements(12);
    // ta.particles.X(0)=VECTOR<double,2>(6,1);
    // ta.particles.X(1)=VECTOR<double,2>(8,2);
    // ta.particles.X(2)=VECTOR<double,2>(4,2);
    // ta.particles.X(3)=VECTOR<double,2>(8,4);
    // ta.particles.X(4)=VECTOR<double,2>(4,4);
    // ta.particles.X(5)=VECTOR<double,2>(1,3);
    // ta.particles.X(6)=VECTOR<double,2>(10,5);
    // ta.particles.X(7)=VECTOR<double,2>(9,9);
    // ta.particles.X(8)=VECTOR<double,2>(7,7);
    // ta.particles.X(9)=VECTOR<double,2>(5,10);
    // ta.particles.X(10)=VECTOR<double,2>(3,8);
    // ta.particles.X(11)=VECTOR<double,2>(2,6);
    // ta.mesh.number_nodes=12;
    // ta.mesh.elements.Exact_Resize(13);
    // ta.mesh.elements(0).Set(4,3,8);
    // ta.mesh.elements(1).Set(3,7,8);
    // ta.mesh.elements(2).Set(1,6,3);
    // ta.mesh.elements(3).Set(2,1,3);
    // ta.mesh.elements(4).Set(0,1,2);
    // ta.mesh.elements(5).Set(4,8,10);
    // ta.mesh.elements(6).Set(10,8,9);
    // ta.mesh.elements(7).Set(8,7,9);
    // ta.mesh.elements(8).Set(5,4,2);
    // ta.mesh.elements(9).Set(11,10,4);
    // ta.mesh.elements(10).Set(5,4,11);
    // ta.mesh.elements(11).Set(2,4,3);
    // ta.mesh.elements(12).Set(3,7,6);
    // VORONOI_2D<double> vr;
    // vr.Initialize_With_A_Triangulated_Area(ta);
    // exit(0);

    // ./mpm_2d -test 4 -gres 32 -pn 400000 -exclude 0.01 -delaunay -delaunay_maxl 0.03 -delaunay_mina 10 -voronoi -dt 1e-4 -fj 50


    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse(true);
    Run_Simulation<VECTOR<double,2> >(parse_args);
    LOG::Finish_Logging();
    return 0;
}


