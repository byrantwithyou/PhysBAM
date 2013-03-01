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
#include "VORONOI_2D.h"
#include "MPM_SIMULATION.h"

using namespace PhysBAM;

template<class T>
void Initialize(int test,MPM_SIMULATION<VECTOR<T,2> >& sim,VORONOI_2D<T>& voronoi,PARSE_ARGS& parse_args,int grid_res,int particle_res,int particle_count,T density_scale,T ym,T pr)
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    parse_args.Parse();

    T object_mass=1;

    // geometry setting
    switch(test){
        case 1: // stretching beam
            sim.grid.Initialize(TV_INT(2*grid_res+1,0.5*grid_res+1),RANGE<TV>(TV(-1,-0.25),TV(1,0.25)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.4*particle_res+1,0.24*particle_res+1),RANGE<TV>(TV(-0.2,-0.12),TV(0.2,0.12)));
            sim.ground_level=-100;
            sim.dirichlet_box.Append(RANGE<TV>(TV(0.16,-0.15),TV(0.25,0.15)));
            sim.dirichlet_velocity.Append(TV(0.2,0));
            sim.dirichlet_box.Append(RANGE<TV>(TV(-0.25,-0.15),TV(-0.16,0.15)));
            sim.dirichlet_velocity.Append(TV(-0.2,0));
            voronoi.Initialize_With_A_Regular_Grid_Of_Particles(GRID<TV>(TV_INT(0.4*particle_res+1,0.24*particle_res+1),RANGE<TV>(TV(-0.2,-0.12),TV(0.2,0.12))));
            break;
        case 2: // shit fall
            sim.grid.Initialize(TV_INT(1*grid_res+1,1*grid_res+1),RANGE<TV>(TV(-0.5,-0.7),TV(0.5,0.3)));
            sim.particles.Initialize_X_As_A_Grid(TV_INT(0.2*particle_res+1,0.2*particle_res+1),RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1)));
            sim.ground_level=-0.3;
            // voronoi.Initialize_With_A_Regular_Grid_Of_Particles(GRID<TV>(TV_INT(0.2*particle_res+1,0.2*particle_res+1),RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1))));
            break;
        default: PHYSBAM_FATAL_ERROR("Missing test");};

    // material setting
    sim.Initialize();
    if(abs(sim.grid.dX(0)-sim.grid.dX(1))>(T)1e-10){
        LOG::cout<<"grid dx: "<<sim.grid.dX<<std::endl;
        exit(0);}
    switch(test){
        case 1:
            object_mass=(T)1200*density_scale*RANGE<TV>(TV(-0.2,-0.12),TV(0.2,0.12)).Size();
            ym*=(T)1e3;
            for(int p=0;p<sim.particles.number;p++){
                T this_ym=0.7*ym/0.04*sqr(sim.particles.X(p)(0))+0.3*ym;
                sim.mu(p)=(this_ym/((T)2*((T)1+pr)));
                sim.lambda(p)=(this_ym*pr/(((T)1+pr)*((T)1-2*pr)));}
            sim.use_gravity=false;
            break;
        case 2:
            object_mass=(T)1200*density_scale*(RANGE<TV>(TV(-0.1,-0.1),TV(0.1,0.1))).Size();
            ym*=(T)1e5;
            sim.mu.Fill(ym/((T)2*((T)1+pr)));
            sim.lambda.Fill(ym*pr/(((T)1+pr)*((T)1-2*pr)));
            sim.use_visco_plasticity=true;
            sim.visco_nu=1e4;
            sim.visco_tau=1000;
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
    VORONOI_2D<T> voronoi;

    int test_number=-1;
    std::string output_directory="";
    bool use_output_directory=false;
    sim.dt=(T)1e-3;
    int frame_jump=20;
    int grid_res=16,particle_res=32,particle_count=100;
    T density_scale=(T)1;
    T ym_scale=1;
    T pr=(T).3;
    sim.xi=(T)0;
    parse_args.Add("-test",&test_number,"test","test number");
    parse_args.Add("-o",&output_directory,&use_output_directory,"o","output directory");
    parse_args.Add("-dt",&sim.dt,"dt","dt");
    parse_args.Add("-fj",&frame_jump,"fj","frame jump");
    parse_args.Add("-dump_matrix",&sim.dump_matrix,"dump linear system");
    parse_args.Add("-test_system",&sim.test_system,"test linear system");
    parse_args.Add("-stiffness",&ym_scale,"value","scale stiffness");
    parse_args.Add("-poisson_ratio",&pr,"value","poisson's ratio");
    parse_args.Add("-hardening",&sim.xi,"value","harderning coefficient");
    parse_args.Add("-flip",&sim.FLIP_alpha,"value","flip fraction");
    parse_args.Add("-gres",&grid_res,"value","grid resolution");
    parse_args.Add("-pres",&particle_res,"value","particle resolution");
    parse_args.Add("-pn",&particle_count,"value","particle number");
    parse_args.Add("-rho",&density_scale,"value","scale object density");
    parse_args.Parse(true);

    typedef VECTOR<int,TV::m> TV_INT;

    if(!use_output_directory) output_directory=STRING_UTILITIES::string_sprintf("MPM_%dD_test_%d",TV::m,test_number);

    Initialize(test_number,sim,voronoi,parse_args,grid_res,particle_res,particle_count,density_scale,ym_scale,pr);

    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,output_directory);
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
    voronoi.Build_Boundary_Segments();
    for(int i=0;i<voronoi.X.m;i++){
        if(voronoi.type(i)==1) Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,0,0));
        else if(voronoi.type(i)==2)Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(1,1,0));
        else Add_Debug_Particle(voronoi.X(i),VECTOR<T,3>(0,0,1));}
    for(int s=0;s<voronoi.boundary_segments.m;s++){
        int i,j;voronoi.boundary_segments(s).Get(i,j);
        Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.boundary_segments(s))),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));}
    for(int b=0;b<sim.rigid_ball.m;b++){
        for(int k=0;k<50;k++){
            T theta=k*2.0*3.14/50.0;
            TV disp;disp(0)=sim.rigid_ball(b).radius*cos(theta);disp(1)=sim.rigid_ball(b).radius*sin(theta);
            Add_Debug_Particle(sim.rigid_ball(b).center+disp,VECTOR<T,3>(1,0,0));}}
    for(int b=0;b<sim.rigid_box.m;b++){
        Add_Debug_Object(VECTOR<TV,TV::m>(sim.rigid_box(b).corner,sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(0)),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));
        Add_Debug_Object(VECTOR<TV,TV::m>(sim.rigid_box(b).corner,sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(1)),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));
        Add_Debug_Object(VECTOR<TV,TV::m>(sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(0)+sim.rigid_box(b).edges.Column(1),sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(0)),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));
        Add_Debug_Object(VECTOR<TV,TV::m>(sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(0)+sim.rigid_box(b).edges.Column(1),sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(1)),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));}
    for(int i=0;i<50;i++) Add_Debug_Particle(TV(-1+i*0.04,sim.ground_level),VECTOR<T,3>(0,0,1));

    Flush_Frame<TV>("mpm");

    for(int f=1;f<29977;f++){
        if(f==12500) sim.dirichlet_velocity.Fill(TV());
        TIMING_START;
        LOG::cout<<"TIMESTEP "<<f<<std::endl;
        sim.Advance_One_Time_Step_Backward_Euler();
        TIMING_END("Current time step totally");
        if(f%frame_jump==0){
            // for(int i=0;i<sim.particles.X.m;i++) if(sim.valid(i)) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
            // for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),sim.grid.counts));it.Valid();it.Next()){
            //     Add_Debug_Particle(sim.grid.Node(it.index),VECTOR<T,3>(1,1,1));
            //     Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,sim.node_force(it.index));}
            voronoi.Crack(sim.particles.X,sim.grid.dX.Min()*1.3);
            voronoi.Build_Association();
            voronoi.Build_Boundary_Segments();
            voronoi.Deform_Mesh_Using_Particle_Deformation(sim.particles.Xm,sim.particles.X,sim.particles.Fe,sim.particles.Fp);
            for(int s=0;s<voronoi.boundary_segments.m;s++){
                int i,j;voronoi.boundary_segments(s).Get(i,j);
                Add_Debug_Object(VECTOR<TV,TV::m>(voronoi.X.Subset(voronoi.boundary_segments(s))),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));}
            for(int b=0;b<sim.rigid_ball.m;b++){
                for(int k=0;k<50;k++){
                    T theta=k*2.0*3.14/50.0;
                    TV disp;disp(0)=sim.rigid_ball(b).radius*cos(theta);disp(1)=sim.rigid_ball(b).radius*sin(theta);
                    Add_Debug_Particle(sim.rigid_ball(b).center+disp,VECTOR<T,3>(1,0,0));}}
            for(int b=0;b<sim.rigid_box.m;b++){
                Add_Debug_Object(VECTOR<TV,TV::m>(sim.rigid_box(b).corner,sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(0)),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));
                Add_Debug_Object(VECTOR<TV,TV::m>(sim.rigid_box(b).corner,sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(1)),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));
                Add_Debug_Object(VECTOR<TV,TV::m>(sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(0)+sim.rigid_box(b).edges.Column(1),sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(0)),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));
                Add_Debug_Object(VECTOR<TV,TV::m>(sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(0)+sim.rigid_box(b).edges.Column(1),sim.rigid_box(b).corner+sim.rigid_box(b).edges.Column(1)),VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,0,0));}
            for(int i=0;i<50;i++) Add_Debug_Particle(TV(-1+i*0.04,sim.ground_level),VECTOR<T,3>(0,0,1));
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


