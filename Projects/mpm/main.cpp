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
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <omp.h>
#include "TIMING.h"
#include "MPM_SIMULATION.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    static const int dimension=2;
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,dimension> TV;
    typedef VECTOR<int,dimension> TV_INT;

    MPM_SIMULATION<TV> sim;

    static const int test=4;
    if(test==1){ // cube falling on ground
        static const int grid_res=64;
        TV_INT grid_counts(1*grid_res,1*grid_res);
        RANGE<TV> grid_box(TV(-0.5,-0.5),TV(0.5,0.5));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=300;
        TV_INT particle_counts(0.2*particle_res,0.2*particle_res);
        RANGE<TV> particle_box(TV(-0.1,-0.1),TV(0.1,0.1));
        sim.particles.Initialize_X_As_A_Grid(particle_counts,particle_box);
        T object_mass=400*0.04;
        for(int p=0;p<sim.particles.number;p++){
            sim.particles.V(p)=TV();
            sim.particles.mass(p)=object_mass/sim.particles.number;
            sim.particles.Fe(p)=MATRIX<T,TV::m>::Identity_Matrix();
            sim.particles.Fp(p)=MATRIX<T,TV::m>::Identity_Matrix();}
        sim.dt=1e-3;
        T ym=3000;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity=false;
        sim.use_gravity=true;
        sim.ground_level=-0.3;
        sim.theta_c=0.025;
        sim.theta_s=0.0075;
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
        sim.dt=1e-4;
        T ym=1500;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity=false;
        sim.use_gravity=true;
        sim.ground_level=-100;
        sim.theta_c=0.025;
        sim.theta_s=0.0075;
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
        sim.dt=1e-4;
        T ym=500;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity=false;
        sim.use_gravity=false;
        sim.ground_level=-100;
        sim.theta_c=0.025;
        sim.theta_s=0.0075;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}
    if(test==4){ // cube falling on beam
        static const int grid_res=64;
        TV_INT grid_counts(0.7*grid_res,0.85*grid_res);
        RANGE<TV> grid_box(TV(-0.35,-0.3),TV(0.35,0.55));
        GRID<TV> grid(grid_counts,grid_box);
        sim.grid=grid;
        static const int particle_res=200;
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
        sim.dt=1e-4;
        T ym=2500;
        T pr=0.3;
        sim.mu0=ym/((T)2*((T)1+pr));
        sim.lambda0=ym*pr/(((T)1+pr)*((T)1-2*pr));
        sim.xi=0;
        sim.use_plasticity=false;
        sim.use_gravity=true;
        sim.ground_level=-100;
        sim.theta_c=0.025;
        sim.theta_s=0.0075;
        sim.FLIP_alpha=0.95;
        sim.friction_coefficient=0.6;}
    

    sim.Initialize();

    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,"output");
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
    Flush_Frame<TV>("mpm");

    for(int f=1;f<1000000;f++){
        TIMING_START;
        LOG::cout<<"TIMESTEP "<<f<<std::endl;
        sim.Advance_One_Time_Step_Backward_Euler();
        TIMING_END("Current time step totally");
        if(f%100==0){
            for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
            Flush_Frame<TV>("mpm");}
        LOG::cout<<std::endl;
    }

    return 0;
}
