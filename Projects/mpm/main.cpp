//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
// tests
//#####################################################################
//   1. cube falling on ground
//   2. bending beam
//   3. stretching beam
//   4. cube falling on beam
//   5. notch test
//   6. wall test
//   7. dropping sphere
//   8. two ring hit each other
//#####################################################################
/*
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

using namespace PhysBAM;

template<class T>
void Initialize(int test,MPM_SIMULATION<VECTOR<T,2> >& sim,PARSE_ARGS& parse_args,int grid_res,int particle_res,T object_mass)
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    parse_args.Parse();
    sim.ym0*=3000;
    switch(test){
        case 1: // cube falling on ground
            sim.grid.Initialize(TV_INT(1*grid_res,1.5*grid_res),RANGE<TV>(TV(-0.5,-1.0),TV(0.5,0.5)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.2*particle_res,0.2*particle_res),RANGE<TV>::Centered_Box()*(T).1);
            sim.use_plasticity_yield=false;
            sim.ground_level=-0.7;
            break;
        case 2: // bending beam
            sim.grid.Initialize(TV_INT(1*grid_res,1*grid_res),RANGE<TV>(TV(-0.5,-0.5),TV(0.5,0.5)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.6*particle_res,0.2*particle_res),RANGE<TV>(TV(-0.3,-0.1),TV(0.3,0.1)));
            sim.use_plasticity_yield=false;
            sim.ground_level=-100;
            sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(-0.25,10)));
            sim.dirichlet_velocity.Append(TV());
            break;
        case 3: // stretching beam
            sim.grid.Initialize(TV_INT(4*grid_res,0.6*grid_res),RANGE<TV>(TV(-2,-0.3),TV(2,0.3)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.6*particle_res,0.2*particle_res),RANGE<TV>(TV(-0.3,-0.1),TV(0.3,0.1)));
            sim.use_plasticity_yield=false;
            sim.ground_level=-100;
            sim.use_gravity=false;
            sim.dirichlet_box.Append(RANGE<TV>(TV(0.25,-10),TV(10,10)));
            sim.dirichlet_velocity.Append(TV(0.2,0));
            sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(-0.25,10)));
            sim.dirichlet_velocity.Append(TV(-0.2,0));
            break;
        case 4: // cube falling on beam
            sim.grid.Initialize(TV_INT(1.0*grid_res,1.2*grid_res),RANGE<TV>(TV(-0.5,-0.5),TV(0.5,0.7)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.6*particle_res,0.2*particle_res),RANGE<TV>(TV(-0.3,-0.1),TV(0.3,0.1)));
            sim.particles.Add_X_As_A_Grid(TV_INT(0.2*particle_res,0.2*particle_res),RANGE<TV>(TV(-0.1,0.3),TV(0.1,0.5)));
            sim.use_plasticity_yield=false;
            sim.ground_level=-100;
            sim.dirichlet_box.Append(RANGE<TV>(TV(0.25,-10),TV(10,10)));
            sim.dirichlet_velocity.Append(TV());
            sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(-0.25,10)));
            sim.dirichlet_velocity.Append(TV());
            break;
        case 5: // notch test
            sim.grid.Initialize(TV_INT(0.4*grid_res,0.8*grid_res),RANGE<TV>(TV(-0.2,-0.4),TV(0.2,0.4)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.2*particle_res,0.3*particle_res),RANGE<TV>(TV(-0.1,-0.15),TV(0.1,0.15)));
            sim.particles.Reduce_X_As_A_Ball(RANGE<TV>(TV(0.07,-0.03),TV(0.13,0.03)));
            sim.ground_level=-100;
            sim.dirichlet_box.Append(RANGE<TV>(TV(-10,0.13),TV(10,10)));
            sim.dirichlet_velocity.Append(TV(0,0.2));
            sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(10,-0.13)));
            sim.dirichlet_velocity.Append(TV(0,-0.2));
            sim.yield_max=1.1;
            sim.yield_min=1.0/sim.yield_max;
            sim.use_plasticity_clamp=true;
            sim.clamp_max=1.3;
            sim.clamp_min=1.0/sim.clamp_max;
            break;
        case 6: // wall test
            sim.grid.Initialize(TV_INT(3.2*grid_res,0.95*grid_res),RANGE<TV>(TV(-0.2,-0.45),TV(3.0,0.5)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.1*particle_res,0.8*particle_res),RANGE<TV>(TV(-0.05,-0.4),TV(0.05,0.4)));
            sim.ground_level=-0.4;
            // sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10),TV(10,-0.38)));
            // sim.dirichlet_velocity.Append(TV());
            sim.rigid_ball.Append(SPHERE<TV>(TV(-0.25,0),0.03));
            sim.rigid_ball_velocity.Append(TV(1,0));
            sim.yield_max=1.1;
            sim.yield_min=1.0/sim.yield_max;
            sim.use_plasticity_clamp=true;
            sim.clamp_max=1.3;
            sim.clamp_min=1.0/sim.clamp_max;
            break;
        case 7: // dropping sphere
            sim.grid.Initialize(TV_INT(0.4*grid_res,0.45*grid_res),RANGE<TV>(TV(-0.2,-0.25),TV(0.2,0.2)));
            sim.particles.Initialize_X_As_A_Ball(TV_INT(0.2*particle_res,0.2*particle_res),RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1)));
            sim.particles.Reduce_X_As_A_Ball(RANGE<TV>(TV(-0.08,-0.08),TV(0.08,0.08)));
            sim.use_plasticity_yield=false;
            sim.ground_level=-0.2;
            sim.yield_max=1.1;
            sim.yield_min=1.0/sim.yield_max;
            sim.use_plasticity_clamp=false;
            sim.clamp_max=1.3;
            sim.clamp_min=1.0/sim.clamp_max;
            break;
        case 8: // two ring hitting each other
            sim.grid.Initialize(TV_INT(1.0*grid_res,0.4*grid_res),RANGE<TV>(TV(-0.5,-0.2),TV(0.5,0.2)));
            sim.particles.Initialize_X_As_A_Ball(TV_INT(0.2*particle_res,0.2*particle_res),RANGE<TV>(TV(-0.3,-0.1),TV(-0.1,0.1)));
            sim.particles.Reduce_X_As_A_Ball(RANGE<TV>(TV(-0.28,-0.08),TV(-0.12,0.08)));
            sim.particles.Add_X_As_A_Ball(TV_INT(0.2*particle_res,0.2*particle_res),RANGE<TV>(TV(0.1,-0.1),TV(0.3,0.1)));
            sim.ground_level=-100;
            sim.particles.Reduce_X_As_A_Ball(RANGE<TV>(TV(0.12,-0.08),TV(0.28,0.08)));
            sim.use_plasticity_yield=false;
            sim.yield_max=1.1;
            sim.yield_min=1.0/sim.yield_max;
            sim.use_plasticity_clamp=false;
            sim.clamp_max=1.3;
            sim.clamp_min=1.0/sim.clamp_max;
            sim.use_gravity=false;
            for(int p=0;p<sim.particles.number;p++){
                if(sim.particles.X(p)(0)<0) sim.particles.V(p)=TV(2,0);
                else sim.particles.V(p)=TV(-2,0);}
            break;
        default: PHYSBAM_FATAL_ERROR("Missing test");};

    for(int p=0;p<sim.particles.number;p++){
        sim.particles.mass(p)=object_mass/sim.particles.number;
        sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
        sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
}

template<class T>
void Initialize(int test,MPM_SIMULATION<VECTOR<T,3> >& sim,PARSE_ARGS& parse_args,int grid_res,int particle_res,T object_mass)
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    parse_args.Parse();
    sim.ym0*=3000;
    switch(test){
        case 4: // cube falling on beam
            sim.grid.Initialize(TV_INT(0.7*grid_res,0.85*grid_res,0.4*grid_res),RANGE<TV>(TV(-0.35,-0.3,-0.2),TV(0.35,0.55,0.2)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.6*particle_res,0.2*particle_res,0.2*particle_res),RANGE<TV>(TV(-0.3,-0.1,-0.1),TV(0.3,0.1,0.1)));
            sim.particles.Add_X_As_A_Grid(TV_INT(0.2*particle_res,0.2*particle_res,0.2*particle_res),RANGE<TV>(TV(-0.1,0.3,-0.1),TV(0.1,0.5,0.1)));
            sim.dirichlet_box.Append(RANGE<TV>(TV(0.25,-10,-10),TV(10,10,10)));
            sim.dirichlet_velocity.Append(TV());
            sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10,-10),TV(-0.25,10,10)));
            sim.dirichlet_velocity.Append(TV());
            sim.use_plasticity_yield=false;
            sim.ground_level=-100;
            break;
        case 6: // wall test
            sim.grid.Initialize(TV_INT(0.4*grid_res,0.95*grid_res,0.95*grid_res),RANGE<TV>(TV(-0.2,-0.45,-0.45),TV(0.2,0.5,0.5)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.1*particle_res,0.8*particle_res,0.8*particle_res),RANGE<TV>(TV(-0.05,-0.4,-0.4),TV(0.05,0.4,0.4)));
            sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10,-10),TV(10,-0.38,10)));
            sim.dirichlet_velocity.Append(TV());
            sim.rigid_ball.Append(SPHERE<TV>(TV(-0.25,0,0),0.06));
            sim.rigid_ball_velocity.Append(TV(5,0,0));
            sim.yield_max=1.1;
            sim.yield_min=1.0/sim.yield_max;
            sim.use_plasticity_clamp=true;
            sim.clamp_max=1.3;
            sim.clamp_min=1.0/sim.clamp_max;
            sim.ground_level=-0.4;
            break;
        default: PHYSBAM_FATAL_ERROR("Missing test");};

    for(int p=0;p<sim.particles.number;p++){
        sim.particles.mass(p)=object_mass/sim.particles.number;
        sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
        sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
}

template<class TV>
void Run_Simulation(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef float RW;

    MPM_SIMULATION<TV> sim;

    int test_number=-1;
    std::string output_directory="";
    bool use_output_directory=false;
    sim.dt=(T)1e-3;
    int frame_jump=20;
    int grid_res=16,particle_res=32;
    T mass=(T)1;
    sim.ym0=1;
    sim.pr0=(T).3;
    sim.xi=(T)0;
    parse_args.Add("-test",&test_number,"test","test number");
    parse_args.Add("-o",&output_directory,&use_output_directory,"o","output directory");
    parse_args.Add("-dt",&sim.dt,"dt","dt");
    parse_args.Add("-fj",&frame_jump,"fj","frame jump");
    parse_args.Add("-dump_matrix",&sim.dump_matrix,"dump linear system");
    parse_args.Add("-test_system",&sim.test_system,"test linear system");
    parse_args.Add("-stiffness",&sim.ym0,"value","scale stiffness");
    parse_args.Add("-poisson_ratio",&sim.pr0,"value","poisson's ratio");
    parse_args.Add("-hardening",&sim.xi,"value","hardening coefficient");
    parse_args.Add("-flip",&sim.FLIP_alpha,"value","flip fraction");
    parse_args.Add("-gres",&grid_res,"value","grid resolution");
    parse_args.Add("-pres",&particle_res,"value","particle resolution");
    parse_args.Add("-m",&mass,"value","object total mass");
    parse_args.Parse(true);

    typedef VECTOR<int,TV::m> TV_INT;

    if(!use_output_directory) output_directory=STRING_UTILITIES::string_sprintf("MPM_%dD_test_%d",TV::m,test_number);

    Initialize(test_number,sim,parse_args,grid_res,particle_res,mass);

    sim.Initialize();

    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,output_directory);
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
    Flush_Frame<TV>("mpm");

    for(int f=1;f<1000000;f++){
        TIMING_START;
        LOG::cout<<"TIMESTEP "<<f<<std::endl;
        sim.Advance_One_Time_Step_Backward_Euler();
        TIMING_END("Current time step totally");
        if(f%frame_jump==0){
            for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
            if(sim.ground_level>-10){
                

            }


            Flush_Frame<TV>("mpm");}
        LOG::cout<<std::endl;}
}

int main(int argc,char *argv[])
{
    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"Use 3D");
    parse_args.Parse(true);

    if(use_3d)
        Run_Simulation<VECTOR<double,3> >(parse_args);
    else
        Run_Simulation<VECTOR<double,2> >(parse_args);
    return 0;
}
*/
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
#include "SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL.h"

using namespace PhysBAM;
int main(int argc,char *argv[])
{
    typedef double T;
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;

    GEOMETRY_PARTICLES<TV> particles;
    VECTOR<int,TV::m> count(50,50,50);
    RANGE<TV> box(TV(-0.5,-0.5,-0.5),TV(0.5,0.5,0.5));
    GRID<TV> grid(count,box);
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+count));it.Valid();it.Next()){
        TV x=grid.X(it.index);
        sample_X.Append(x);}
    particles.Resize(sample_X.m);
    particles.X=sample_X;
    LOG::cout<<"wtf"<<std::endl;
    SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL<TV> sr;
    ARRAY<TV> Xbar;
    ARRAY<MATRIX<T,TV::m> > G;
    sr.Compute_Kernal_Centers_And_Transformation(particles,0.03,0.06,0.95,15,4,1400,0.5,Xbar,G);

    return 0;
}
