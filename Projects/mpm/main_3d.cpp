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
#include "MPM_SURFACE_2D.h"
#include "MPM_SIMULATION.h"

using namespace PhysBAM;

template<class T>
void Initialize(int test,MPM_SIMULATION<VECTOR<T,3> >& sim,PARSE_ARGS& parse_args,int grid_res,int particle_res,int particle_count,T object_mass)
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
        case 5: // notch test
            sim.grid.Initialize(TV_INT(0.4*grid_res,0.8*grid_res,0.4*grid_res),RANGE<TV>(TV(-0.2,-0.4,-0.2),TV(0.2,0.4,0.2)));
            sim.particles.Initialize_X_As_A_Randomly_Sampled_Box(particle_count,RANGE<TV>(TV(-0.1,-0.2,-0.05),TV(0.1,0.2,0.05)));
            {CYLINDER<T> cylinder(TV(0.1,0,-10),TV(0.1,0,10),0.03);
                ARRAY<TV> old_X;
                old_X.Resize(sim.particles.X.m);
                for(int i=0;i<old_X.m;i++) old_X(i)=sim.particles.X(i);
                ARRAY<TV> sample_X;
                for(int i=0;i<old_X.m;i++)
                    if(!cylinder.Lazy_Inside(old_X(i)))
                        sample_X.Append(old_X(i));
                sim.particles.Resize(sample_X.m);
                sim.particles.X=sample_X;}
            sim.ground_level=-100;
            sim.dirichlet_box.Append(RANGE<TV>(TV(-10,0.19,-10),TV(10,10,10)));
            sim.dirichlet_velocity.Append(TV(0,0.1,0));
            sim.dirichlet_box.Append(RANGE<TV>(TV(-10,-10,-10),TV(10,-0.19,10)));
            sim.dirichlet_velocity.Append(TV(0,-0.1,0));
            sim.use_gravity=false;
            sim.use_plasticity_yield=true;
            sim.yield_max=2.0;
            sim.yield_min=-100;
            sim.use_plasticity_clamp=false;
            sim.clamp_max=1.8;
            sim.clamp_min=1.0/sim.clamp_max;
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
    int grid_res=16,particle_res=32,particle_count=100;
    T mass=(T)40;
    sim.ym0=1;
    sim.pr0=(T).3;
    parse_args.Add("-test",&test_number,"test","test number");
    parse_args.Add("-o",&output_directory,&use_output_directory,"o","output directory");
    parse_args.Add("-dt",&sim.dt,"dt","dt");
    parse_args.Add("-fj",&frame_jump,"fj","frame jump");
    parse_args.Add("-dump_matrix",&sim.dump_matrix,"dump linear system");
    parse_args.Add("-test_system",&sim.test_system,"test linear system");
    parse_args.Add("-stiffness",&sim.ym0,"value","scale stiffness");
    parse_args.Add("-poisson_ratio",&sim.pr0,"value","poisson's ratio");
    parse_args.Add("-flip",&sim.FLIP_alpha,"value","flip fraction");
    parse_args.Add("-gres",&grid_res,"value","grid resolution");
    parse_args.Add("-pres",&particle_res,"value","particle resolution");
    parse_args.Add("-pn",&particle_count,"value","particle number");
    parse_args.Add("-m",&mass,"value","object total mass");
    parse_args.Parse(true);

    typedef VECTOR<int,TV::m> TV_INT;

    if(!use_output_directory) output_directory=STRING_UTILITIES::string_sprintf("MPM_%dD_test_%d",TV::m,test_number);

    Initialize(test_number,sim,parse_args,grid_res,particle_res,particle_count,mass);

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
            for(int i=0;i<sim.particles.X.m;i++) if(sim.valid(i)) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
            Flush_Frame<TV>("mpm");}
        LOG::cout<<std::endl;}
}

int main(int argc,char *argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse(true);
    Run_Simulation<VECTOR<double,3> >(parse_args);
    return 0;
}
