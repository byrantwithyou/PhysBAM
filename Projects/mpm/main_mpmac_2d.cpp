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
#include "MPMAC.h"
#include "TIMING.h"
#include <omp.h>
using namespace PhysBAM;

template<class TV>
void Run_Simulation(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,2> TV_INT;
    typedef float RW;
    MPMAC<TV> sim;
    int test_number=-1;
    std::string output_directory="";
    bool use_output_directory=false;
    sim.dt=(T)1e-3;
    sim.FLIP_alpha=0.95;
    int frame_jump=20;
    int grid_res=16,particle_res=32,particle_count=100;
    T particle_exclude_radius=0;
    
    parse_args.Add("-test",&test_number,"test","test number");
    parse_args.Add("-o",&output_directory,&use_output_directory,"o","output directory");
    parse_args.Add("-dt",&sim.dt,"dt","dt");
    parse_args.Add("-fj",&frame_jump,"fj","frame jump");
    parse_args.Add("-flip",&sim.FLIP_alpha,"value","flip fraction");
    parse_args.Add("-gres",&grid_res,"value","grid resolution");
    parse_args.Add("-pres",&particle_res,"value","particle resolution");
    parse_args.Add("-pn",&particle_count,"value","particle number");
    parse_args.Add("-exclude",&particle_exclude_radius,"value","particle exclude radius when using pn");
    parse_args.Parse(true);
    parse_args.Parse();

    if(!use_output_directory) output_directory=STRING_UTILITIES::string_sprintf("MPMAC_%dD_test_%d",TV::m,test_number);

    // geometry setting
    switch(test_number){
        case 1:{ // all materials together
            sim.grid.Initialize(TV_INT(10.2*grid_res+1,3.2*grid_res+1),RANGE<TV>(TV(-1.6,-1.6),TV(8.6,1.6)));
            sim.uniform_density=false;
            
            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-1.48,-1.48),TV(8.48,-1.3)),particle_exclude_radius);
            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-1.48,-1.3),TV(-0.5,0.4)),particle_exclude_radius);
            int c1=sim.particles.number;
            sim.particles.Set_Material_Properties(0,c1,
                (T)8/1000, // mass per particle
                0, // mu
                0, // lambda
                false); // compress
            sim.particles.Set_Plasticity(0,c1,
                false,-1000,1.2, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(0,c1,
                false,100, // visco_nu
                5000, // visco_tau
                0); // visco_kappa
            sim.particles.Set_Initial_State(0,c1,
                MATRIX<T,TV::m>::Identity_Matrix(), // Fe
                MATRIX<T,TV::m>::Identity_Matrix(), // Fp
                TV(0,0)); // initial velocity
            
//            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-0.5,-1.45),TV(0.5,0)),particle_exclude_radius);
//            int c2=sim.particles.number-c1;
//            sim.particles.Set_Material_Properties(c1,c2,
//                                                  (T)10/1000, // mass per particle
//                                                  0, // mu
//                                                  0, // lambda
//                                                  false); // compress
//            sim.particles.Set_Plasticity(c1,c2,
//                                         false,-1000,1.2, // plasticity_yield
//                                         false,-1,1); // plasticity_clamp
//            sim.particles.Set_Visco_Plasticity(c1,c2,
//                                               false,100, // visco_nu
//                                               5000, // visco_tau
//                                               0); // visco_kappa
//            sim.particles.Set_Initial_State(c1,c2,
//                                            MATRIX<T,TV::m>::Identity_Matrix(), // Fe
//                                            MATRIX<T,TV::m>::Identity_Matrix(), // Fp
//                                            TV(0,0)); // initial velocity
//            
//            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(0.8,-1.45),TV(1.45,0)),particle_exclude_radius);
//            int c3=sim.particles.number-c1-c2;
//            sim.particles.Set_Material_Properties(c1+c2,c3,
//                                                  (T)18/1000, // mass per particle
//                                                  0, // mu
//                                                  0, // lambda
//                                                  false); // compress
//            sim.particles.Set_Plasticity(c1+c2,c3,
//                                         false,-1000,1.2, // plasticity_yield
//                                         false,-1,1); // plasticity_clamp
//            sim.particles.Set_Visco_Plasticity(c1+c2,c3,
//                                               false,100, // visco_nu
//                                               5000, // visco_tau
//                                               0); // visco_kappa
//            sim.particles.Set_Initial_State(c1+c2,c3,
//                                            MATRIX<T,TV::m>::Identity_Matrix(), // Fe
//                                            MATRIX<T,TV::m>::Identity_Matrix(), // Fp
//                                            TV(0,0)); // initial velocity
            
            sim.use_gravity=true;

            break;}
     
        default: PHYSBAM_FATAL_ERROR("Missing test");};

    sim.Initialize();
    if(abs(sim.grid.dX(0)-sim.grid.dX(1))>(T)1e-10){
        LOG::cout<<"grid not uniform! dx: "<<sim.grid.dX<<std::endl;
        exit(0);}

    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,output_directory);

    // draw MPM particles
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,1,1));

    Flush_Frame<TV>("mpmac");

    for(int f=1;f<2900977;f++){
        TIMING_START;
        LOG::cout<<"MPM TIMESTEP "<<f<<std::endl;
        sim.Reinitialize();
        sim.Weights();
        LOG::cout<<"Total momentum on particles before rasterize: "<<sim.Get_Total_Momentum_On_Particles()<<std::endl;
        sim.Rasterize();
        LOG::cout<<"Total momentum on faces after rasterize: "<<sim.Get_Total_Momentum_On_Faces()<<std::endl;
        sim.Advection();
        LOG::cout<<"Total momentum on faces after advection: "<<sim.Get_Total_Momentum_On_Faces()<<std::endl;
        sim.Identify_Dirichlet();
        sim.Identify_Neumann();
        for(int x=3;x<sim.mac_grid.counts.x-3;x++){
            sim.cell_neumann(TV_INT(x,3))=true;
            sim.cell_dirichlet(TV_INT(x,3))=false;
            sim.neumann_cell_normal_axis(TV_INT(x,3))=2;}
        for(int x=3;x<sim.mac_grid.counts.x-3;x++){
            sim.cell_neumann(TV_INT(x,sim.mac_grid.counts.y-4))=true;
            sim.cell_dirichlet(TV_INT(x,sim.mac_grid.counts.y-4))=false;
            sim.neumann_cell_normal_axis(TV_INT(x,sim.mac_grid.counts.y-4))=-2;}
        for(int y=3;y<sim.mac_grid.counts.y-3;y++){
            sim.cell_neumann(TV_INT(3,y))=true;
            sim.cell_dirichlet(TV_INT(3,y))=false;
            sim.neumann_cell_normal_axis(TV_INT(3,y))=1;}
        for(int y=3;y<sim.mac_grid.counts.y-3;y++){
            sim.cell_neumann(TV_INT(sim.mac_grid.counts.x-4,y))=true;
            sim.cell_dirichlet(TV_INT(sim.mac_grid.counts.x-4,y))=false;
            sim.neumann_cell_normal_axis(TV_INT(sim.mac_grid.counts.x-4,y))=-1;}
        
//        for(int x=3;x<sim.mac_grid.counts.x/2;x++){
//            for(int y=sim.mac_grid.counts.y/2-3;y<=sim.mac_grid.counts.y/2+3;y++){
//                sim.cell_neumann(TV_INT(x,y))=true;
//                sim.cell_dirichlet(TV_INT(x,y))=false;
//                sim.neumann_cell_normal_axis(TV_INT(x,y))=2;}}
//        
//        for(int x=sim.mac_grid.counts.x*2/3;x<=sim.mac_grid.counts.x-4;x++){
//            for(int y=sim.mac_grid.counts.y/2-3;y<=sim.mac_grid.counts.y/2+3;y++){
//                sim.cell_neumann(TV_INT(x,y))=true;
//                sim.cell_dirichlet(TV_INT(x,y))=false;
//                sim.neumann_cell_normal_axis(TV_INT(x,y))=2;}}
        
        
        sim.Build_Velocity_Divergence();
        sim.Solve_For_Pressure();
        sim.Do_Projection();
        LOG::cout<<"Total momentum on faces after projection: "<<sim.Get_Total_Momentum_On_Faces()<<std::endl;
        sim.Update_Particle_Velocities();
        LOG::cout<<"Total momentum on particles after interpolating back: "<<sim.Get_Total_Momentum_On_Particles()<<std::endl;
        sim.Particle_Based_Body_Collisions();
        sim.Update_Particle_Positions();

        sim.frame++;

        TIMING_END("Current MPM time step totally");

        if(f%frame_jump==0){
            for(int i=0;i<sim.particles.X.m;i++){
//                if(sim.particles.Xm(i).x<-0.51)
//                    Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,1,1));
//                else if(sim.particles.Xm(i).x<0.51)
                    Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,1));
//                else
//                    Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,0,1));
            }

            // visualize Neumann cells and Dirichlet cells and Compressible cells
            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),sim.mac_grid.counts));it.Valid();it.Next()){
                // if(sim.cell_dirichlet(it.index))
                //    Add_Debug_Particle(sim.mac_grid.Center(it.index),VECTOR<T,3>(0.3,0.3,0.3));
                if(sim.cell_neumann(it.index)){
                    PHYSBAM_ASSERT(!sim.cell_dirichlet(it.index));
                    Add_Debug_Particle(sim.mac_grid.Center(it.index),VECTOR<T,3>(0,1,0));}}

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


