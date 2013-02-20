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
#include <omp.h>
#include "TIMING.h"
#include "MPM_SIMULATION.h"

#define MPM_3D

#define CHECK_ARG(x,arg,argdefault) if(arg!=argdefault) x=arg;
using namespace PhysBAM;
int main(int argc,char *argv[])
{
#ifdef MPM_2D
    static const int dimension=2;
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,dimension> TV;
    typedef VECTOR<int,dimension> TV_INT;

    PARSE_ARGS parse_args(argc,argv);
    int test_input=-1;;
    std::string output_directory_input("");
    T dt_input=0;
    int frame_jump_input=-1;
    parse_args.Add("-test",&test_input,"test","test number");
    parse_args.Add("-o",&output_directory_input,"o","output directory");
    parse_args.Add("-dt",&dt_input,"dt","dt");
    parse_args.Add("-fj",&frame_jump_input,"fj","frame jump");
    parse_args.Parse(true);

    int test=5;CHECK_ARG(test,test_input,-1);
    std::string output_directory=std::string("MPM_2D_test")+FILE_UTILITIES::Number_To_String(test);CHECK_ARG(output_directory,output_directory_input,"");

    MPM_SIMULATION<TV> sim;
    if(test==1){ // cube falling on ground
        static const int grid_res=256;
        TV_INT grid_counts(1*grid_res,1*grid_res);
        RANGE<TV> grid_box(TV(-0.5,-0.5),TV(0.5,0.5));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=1800;
        TV_INT particle_counts(0.2*particle_res,0.2*particle_res);
        RANGE<TV> particle_box(TV(-0.1,-0.1),TV(0.1,0.1));
        sim.particles.Initialize_X_As_A_Grid(particle_counts,particle_box);
        T object_mass=400*0.04;
        for(int p=0;p<sim.particles.number;p++){
            sim.particles.V(p)=TV();
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dt=1e-4;CHECK_ARG(sim.dt,dt_input,0);
        T ym=3000;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity_yield=false;
        sim.use_gravity=true;
        sim.ground_level=-0.3;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}

    if(test==2){ // bending beam
        static const int grid_res=64;
        TV_INT grid_counts(1*grid_res,1*grid_res);
        RANGE<TV> grid_box(TV(-0.5,-0.5),TV(0.5,0.5));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=200;
        TV_INT particle_counts(0.6*particle_res,0.2*particle_res);
        RANGE<TV> particle_box(TV(-0.3,-0.1),TV(0.3,0.1));
        sim.particles.Initialize_X_As_A_Grid(particle_counts,particle_box);
        T object_mass=40*0.12;
        for(int p=0;p<sim.particles.number;p++){
            sim.particles.V(p)=TV();
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dirichlet_box.Append(RANGE<TV>(TV(0.25,-10),TV(10,10)));
        sim.dirichlet_velocity.Append(TV());
        sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(-0.25,10)));
        sim.dirichlet_velocity.Append(TV());
        sim.dt=1e-4;CHECK_ARG(sim.dt,dt_input,0);
        T ym=1500;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity_yield=false;
        sim.use_gravity=true;
        sim.ground_level=-100;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}

    if(test==3){ // stretching beam
        static const int grid_res=64;
        TV_INT grid_counts(4*grid_res,0.6*grid_res);
        RANGE<TV> grid_box(TV(-2,-0.3),TV(2,0.3));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=200;
        TV_INT particle_counts(0.6*particle_res,0.2*particle_res);
        RANGE<TV> particle_box(TV(-0.3,-0.1),TV(0.3,0.1));
        sim.particles.Initialize_X_As_A_Grid(particle_counts,particle_box);
        T object_mass=40*0.12;
        for(int p=0;p<sim.particles.number;p++){
            sim.particles.V(p)=TV();
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dirichlet_box.Append(RANGE<TV>(TV(0.25,-10),TV(10,10)));
        sim.dirichlet_velocity.Append(TV(0.2,0));
        sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(-0.25,10)));
        sim.dirichlet_velocity.Append(TV(-0.2,0));
        sim.dt=1e-4;CHECK_ARG(sim.dt,dt_input,0);
        T ym=500;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity_yield=false;
        sim.use_gravity=false;
        sim.ground_level=-100;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}

    if(test==4){ // cube falling on beam
        static const int grid_res=256;
        TV_INT grid_counts(0.7*grid_res,0.85*grid_res);
        RANGE<TV> grid_box(TV(-0.35,-0.3),TV(0.35,0.55));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=800;
        TV_INT particle_counts_beam(0.6*particle_res,0.2*particle_res);
        RANGE<TV> particle_box_beam(TV(-0.3,-0.1),TV(0.3,0.1));
        TV_INT particle_counts_cube(0.2*particle_res,0.2*particle_res);
        RANGE<TV> particle_box_cube(TV(-0.1,0.3),TV(0.1,0.5));
        sim.particles.Initialize_X_As_A_Grid(particle_counts_beam,particle_box_beam);
        sim.particles.Add_X_As_A_Grid(particle_counts_cube,particle_box_cube);
        T object_mass=40*(0.12+0.04);
        for(int p=0;p<sim.particles.number;p++){
            sim.particles.V(p)=TV();
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dirichlet_box.Append(RANGE<TV>(TV(0.25,-10),TV(10,10)));
        sim.dirichlet_velocity.Append(TV());
        sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(-0.25,10)));
        sim.dirichlet_velocity.Append(TV());
        sim.dt=1e-4;CHECK_ARG(sim.dt,dt_input,0);
        T ym=2500;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity_yield=false;
        sim.use_gravity=true;
        sim.ground_level=-100;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}

    if(test==5){ // notch test
        static const int grid_res=64;
        TV_INT grid_counts(0.4*grid_res,0.8*grid_res);
        RANGE<TV> grid_box(TV(-0.2,-0.4),TV(0.2,0.4));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=500;
        TV_INT particle_counts_cube(0.2*particle_res,0.3*particle_res);
        RANGE<TV> particle_box_cube(TV(-0.1,-0.15),TV(0.1,0.15));
        RANGE<TV> ball_box(TV(0.07,-0.03),TV(0.13,0.03));
        sim.particles.Initialize_X_As_A_Grid(particle_counts_cube,particle_box_cube);
        sim.particles.Reduce_X_As_A_Ball(ball_box);
        T object_mass=1;;
        for(int p=0;p<sim.particles.number;p++){
            sim.particles.V(p)=TV();
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dirichlet_box.Append(RANGE<TV>(TV(-10,0.13),TV(10,10)));
        sim.dirichlet_velocity.Append(TV(0,0.2));
        sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(10,-0.13)));
        sim.dirichlet_velocity.Append(TV(0,-0.2));
        sim.dt=1e-4;CHECK_ARG(sim.dt,dt_input,0);
        T ym=3000;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity_yield=true;
        sim.yield_max=1.1;
        sim.yield_min=1.0/sim.yield_max;
        sim.use_plasticity_clamp=true;
        sim.clamp_max=1.3;
        sim.clamp_min=1.0/sim.clamp_max;
        sim.use_gravity=true;
        sim.ground_level=-100;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}

    if(test==6){ // wall test
        static const int grid_res=64;
        TV_INT grid_counts(3.2*grid_res,0.95*grid_res);
        RANGE<TV> grid_box(TV(-0.2,-0.45),TV(3.0,0.5));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=500;
        TV_INT particle_counts_cube(0.1*particle_res,0.8*particle_res);
        RANGE<TV> particle_box_cube(TV(-0.05,-0.4),TV(0.05,0.4));
        sim.particles.Initialize_X_As_A_Grid(particle_counts_cube,particle_box_cube);
        T object_mass=10;
        for(int p=0;p<sim.particles.number;p++){
            sim.particles.V(p)=TV();
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(10,-0.38)));
        sim.dirichlet_velocity.Append(TV());
        sim.rigid_ball.Append(SPHERE<TV>(TV(-0.25,0),0.03));
        sim.rigid_ball_velocity.Append(TV(5,0));
        sim.dt=1e-4;CHECK_ARG(sim.dt,dt_input,0);
        T ym=3000;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity_yield=true;
        sim.yield_max=1.1;
        sim.yield_min=1.0/sim.yield_max;
        sim.use_plasticity_clamp=true;
        sim.clamp_max=1.3;
        sim.clamp_min=1.0/sim.clamp_max;
        sim.use_gravity=true;
        sim.ground_level=-0.4;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}

    if(test==7){ // dropping sphere
        static const int grid_res=64;
        TV_INT grid_counts(0.4*grid_res,0.45*grid_res);
        RANGE<TV> grid_box(TV(-0.2,-0.25),TV(0.2,0.2));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=500;
        TV_INT particle_counts_cube(0.2*particle_res,0.2*particle_res);
        RANGE<TV> particle_box_cube(TV(-0.1,-0.1),TV(0.1,0.1));
        sim.particles.Initialize_X_As_A_Ball(particle_counts_cube,particle_box_cube);
        sim.particles.Reduce_X_As_A_Ball(RANGE<TV>(TV(-0.08,-0.08),TV(0.08,0.08)));
        T object_mass=10;
        for(int p=0;p<sim.particles.number;p++){
            sim.particles.V(p)=TV();
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dt=1e-4;CHECK_ARG(sim.dt,dt_input,0);
        T ym=100000;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity_yield=false;
        sim.yield_max=1.1;
        sim.yield_min=1.0/sim.yield_max;
        sim.use_plasticity_clamp=false;
        sim.clamp_max=1.3;
        sim.clamp_min=1.0/sim.clamp_max;
        sim.use_gravity=true;
        sim.ground_level=-0.2;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}

    if(test==8){ // two ring hitting each other
        static const int grid_res=64;
        TV_INT grid_counts(1.0*grid_res,0.4*grid_res);
        RANGE<TV> grid_box(TV(-0.5,-0.2),TV(0.5,0.2));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=500;
        TV_INT particle_counts_ball1(0.2*particle_res,0.2*particle_res);
        RANGE<TV> particle_box_ball1(TV(-0.3,-0.1),TV(-0.1,0.1));
        TV_INT particle_counts_ball2(0.2*particle_res,0.2*particle_res);
        RANGE<TV> particle_box_ball2(TV(0.1,-0.1),TV(0.3,0.1));
        sim.particles.Initialize_X_As_A_Ball(particle_counts_ball1,particle_box_ball1);
        sim.particles.Reduce_X_As_A_Ball(RANGE<TV>(TV(-0.28,-0.08),TV(-0.12,0.08)));
        sim.particles.Add_X_As_A_Ball(particle_counts_ball2,particle_box_ball2);
        sim.particles.Reduce_X_As_A_Ball(RANGE<TV>(TV(0.12,-0.08),TV(0.28,0.08)));
        T object_mass=10;
        for(int p=0;p<sim.particles.number;p++){
            if(sim.particles.X(p)(0)<0) sim.particles.V(p)=TV(2,0);
            else sim.particles.V(p)=TV(-2,0);
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dt=1e-4;CHECK_ARG(sim.dt,dt_input,0);
        T ym=100000;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity_yield=false;
        sim.yield_max=1.1;
        sim.yield_min=1.0/sim.yield_max;
        sim.use_plasticity_clamp=false;
        sim.clamp_max=1.3;
        sim.clamp_min=1.0/sim.clamp_max;
        sim.use_gravity=false;
        sim.ground_level=-100;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}

    sim.Initialize();

    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,output_directory);
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
    Flush_Frame<TV>("mpm");

    int frame_jump=20;CHECK_ARG(frame_jump,frame_jump_input,-1);
    for(int f=1;f<1000000;f++){
        TIMING_START;
        LOG::cout<<"TIMESTEP "<<f<<std::endl;
        sim.Advance_One_Time_Step_Backward_Euler();
        TIMING_END("Current time step totally");
        if(f%frame_jump==0){
            for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
            for(int b=0;b<sim.rigid_ball.m;b++){
                for(T theta=0;theta<2*3.1415;theta+=2*3.1415/20.0){
                    TV p=sim.rigid_ball(b).center;
                    p(0)+=sim.rigid_ball(b).radius*cos(theta);
                    p(1)+=sim.rigid_ball(b).radius*sin(theta);
                    Add_Debug_Particle(p,VECTOR<T,3>(1,0,0));}}
            if(sim.ground_level>-10)
                for(T x=-5;x<5;x+=0.02)
                    Add_Debug_Particle(TV(x,sim.ground_level),VECTOR<T,3>(0,0,1));
            Flush_Frame<TV>("mpm");}
        LOG::cout<<std::endl;
    }
#endif

#ifdef MPM_3D
    static const int dimension=3;
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,dimension> TV;
    typedef VECTOR<int,dimension> TV_INT;

    PARSE_ARGS parse_args(argc,argv);
    int test_input=-1;;
    std::string output_directory_input("");
    T dt_input=0;
    int frame_jump_input=-1;
    parse_args.Add("-test",&test_input,"test","test number");
    parse_args.Add("-o",&output_directory_input,"o","output directory");
    parse_args.Add("-dt",&dt_input,"dt","dt");
    parse_args.Add("-fj",&frame_jump_input,"fj","frame jump");
    parse_args.Parse(true);

    int test=4;CHECK_ARG(test,test_input,-1);
    std::string output_directory=std::string("MPM_3D_test")+FILE_UTILITIES::Number_To_String(test);CHECK_ARG(output_directory,output_directory_input,"");

    MPM_SIMULATION<TV> sim;

    if(test==4){ // cube falling on beam
        static const int grid_res=32;
        TV_INT grid_counts(0.7*grid_res,0.85*grid_res,0.4*grid_res);
        RANGE<TV> grid_box(TV(-0.35,-0.3,-0.2),TV(0.35,0.55,0.2));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=64;
        TV_INT particle_counts_beam(0.6*particle_res,0.2*particle_res,0.2*particle_res);
        RANGE<TV> particle_box_beam(TV(-0.3,-0.1,-0.1),TV(0.3,0.1,0.1));
        TV_INT particle_counts_cube(0.2*particle_res,0.2*particle_res,0.2*particle_res);
        RANGE<TV> particle_box_cube(TV(-0.1,0.3,-0.1),TV(0.1,0.5,0.1));
        sim.particles.Initialize_X_As_A_Grid(particle_counts_beam,particle_box_beam);
        sim.particles.Add_X_As_A_Grid(particle_counts_cube,particle_box_cube);
        T object_mass=1;
        for(int p=0;p<sim.particles.number;p++){
            if(sim.particles.X(p)(1)>0.2) sim.particles.V(p)=TV(0,-1.5,0);
            else sim.particles.V(p)=TV();
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dirichlet_box.Append(RANGE<TV>(TV(0.25,-10,-10),TV(10,10,10)));
        sim.dirichlet_velocity.Append(TV());
        sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10,-10),TV(-0.25,10,10)));
        sim.dirichlet_velocity.Append(TV());
        sim.dt=1e-4;CHECK_ARG(sim.dt,dt_input,0);
        T ym=3000;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity_yield=false;
        sim.use_gravity=true;
        sim.ground_level=-100;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}

    sim.Initialize();

    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,output_directory);
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
    Flush_Frame<TV>("mpm");

    int frame_jump=20;CHECK_ARG(frame_jump,frame_jump_input,-1);
    for(int f=1;f<1000000;f++){
        TIMING_START;
        LOG::cout<<"TIMESTEP "<<f<<std::endl;
        sim.Advance_One_Time_Step_Backward_Euler();
        TIMING_END("Current time step totally");
        if(f%frame_jump==0){
            for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
            Flush_Frame<TV>("mpm");}
        LOG::cout<<std::endl;}

#endif
    return 0;
}
