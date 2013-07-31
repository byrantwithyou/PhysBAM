//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/pow.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Read_Write/READ_WRITE_FORWARD.h>
#include <Tools/Read_Write/TYPED_STREAM.h>
#include <Tools/Utilities/TIMER.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Geometry/Basic_Geometry/POLYGON.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Grids_Uniform_Computations/REINITIALIZATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include "MPMAC.h"
#include "TIMING.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif
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
    sim.uniform_density=false;
    
    parse_args.Add("-test",&test_number,"test","test number");
    parse_args.Add("-o",&output_directory,&use_output_directory,"o","output directory");
    parse_args.Add("-dt",&sim.dt,"dt","dt");
    parse_args.Add("-fj",&frame_jump,"fj","frame jump");
    parse_args.Add("-flip",&sim.FLIP_alpha,"value","flip fraction");
    parse_args.Add("-gres",&grid_res,"value","grid resolution");
    parse_args.Add("-pres",&particle_res,"value","particle resolution");
    parse_args.Add("-pn",&particle_count,"value","particle number");
    parse_args.Add("-exclude",&particle_exclude_radius,"value","particle exclude radius when using pn");
    parse_args.Add("-uniform_rho",&sim.uniform_density,"value","whether use momentum fixed flip");

    parse_args.Parse(true);
    parse_args.Parse();

    if(!use_output_directory) output_directory=STRING_UTILITIES::string_sprintf("MPMAC_%dD_test_%d",TV::m,test_number);

    int air=0;
    // geometry setting
    switch(test_number){
        case 1:{ // dam break
            sim.grid.Initialize(TV_INT(3.2*grid_res+1,3.2*grid_res+1),RANGE<TV>(TV(-1.6,-1.6),TV(1.6,1.6)));
            
            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-1.3,-1.3),TV(-0.8,0.8)),particle_exclude_radius);

            int c1=sim.particles.number;
            sim.particles.Set_Material_Properties(0,c1,
                (T)8/1000, // mass per particle
                0, // mu
                0, // lambda
                false,0); // compress,pressure
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
            
            sim.use_gravity=true;

            break;}
            
        case 2:{ // ralyeigh taylor
            sim.grid.Initialize(TV_INT(3.2*grid_res+1,3.2*grid_res+1),RANGE<TV>(TV(-1.6,-1.6),TV(1.6,1.6)));
            
            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-1.46,-1.46),TV(1.46,-0.5)),particle_exclude_radius);
            int c1=sim.particles.number;
            sim.particles.Set_Material_Properties(0,c1,
                (T)8/1000, // mass per particle
                0, // mu
                0, // lambda
                false,0); // compress,pressure
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
            air=c1;

            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-1.46,-0.49),TV(1.46,0.49)),particle_exclude_radius);
            int c2=sim.particles.number-c1;
            sim.particles.Set_Material_Properties(c1,c2,
                (T)12/1000, // mass per particle
                0, // mu
                0, // lambda
                false,0); // compress,pressure
            sim.particles.Set_Plasticity(c1,c2,
                false,-1000,1.2, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(c1,c2,
                false,100, // visco_nu
                5000, // visco_tau
                0); // visco_kappa
            sim.particles.Set_Initial_State(c1,c2,
                MATRIX<T,TV::m>::Identity_Matrix(), // Fe
                MATRIX<T,TV::m>::Identity_Matrix(), // Fp
                TV(0,0)); // initial velocity
            
            
            sim.use_gravity=true;
            
            break;}
            
        case 3:{ // droplet with air surroundings
            sim.grid.Initialize(TV_INT(3.2*grid_res+1,3.2*grid_res+1),RANGE<TV>(TV(-1.6,-1.6),TV(1.6,1.6)));
            
            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-1.46,-1.46),TV(1.46,1.46)),particle_exclude_radius);
            sim.particles.Delete_Particles_In_Object_Lazy(RANGE<TV>(TV(-1.46,-1.46),TV(1.46,-0.9)));
            sim.particles.Delete_Particles_In_Object_Lazy(SPHERE<TV>(TV(0,1),0.2));
            int c1=sim.particles.number;
            sim.particles.Set_Material_Properties(0,c1,
                (T)1/1000.0, // mass per particle
                0, // mu
                0, // lambda
                false,0); // compress,pressure
            sim.particles.Set_Plasticity(0,c1,
                false,-1000,1.46, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(0,c1,
                false,100, // visco_nu
                5000, // visco_tau
                0); // visco_kappa
            sim.particles.Set_Initial_State(0,c1,
                MATRIX<T,TV::m>::Identity_Matrix(), // Fe
                MATRIX<T,TV::m>::Identity_Matrix(), // Fp
                TV(0,0)); // initial velocity
            
            air=c1;
            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-1.46,-1.46),TV(1.46,-0.9)),particle_exclude_radius);
            sim.particles.Add_Randomly_Sampled_Object(SPHERE<TV>(TV(0,1),0.2),particle_exclude_radius);
            int c2=sim.particles.number-c1;
            sim.particles.Set_Material_Properties(c1,c2,
                (T)1000/1000.0, // mass per particle
                0, // mu
                0, // lambda
                false,0); // compress,pressure
            sim.particles.Set_Plasticity(c1,c2,
                false,-1000,1.2, // plasticity_yield
                false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(c1,c2,
                false,100, // visco_nu
                5000, // visco_tau
                0); // visco_kappa
            sim.particles.Set_Initial_State(c1,c2,
                MATRIX<T,TV::m>::Identity_Matrix(), // Fe
                MATRIX<T,TV::m>::Identity_Matrix(), // Fp
                TV(0,0)); // initial velocity
            
            
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

    T total_time_for_solve=(T)0;
    for(int f=1;f<=26640;f++){
        TIMING_START;
        LOG::cout<<"MPM TIMESTEP "<<f<<std::endl;
        sim.Reinitialize();
        T Timing_reinit=TIMING_GET;
        sim.Weights();
        T Timing_weights=TIMING_GET;
        // LOG::cout<<"Total momentum on particles before rasterize: "<<sim.Get_Total_Momentum_On_Particles()<<std::endl;
        sim.Rasterize();
        T Timing_rasterize=TIMING_GET;
        // LOG::cout<<"Total momentum on faces after rasterize: "<<sim.Get_Total_Momentum_On_Faces()<<std::endl;
        sim.Advection();
        T Timing_advection=TIMING_GET;
        // LOG::cout<<"Total momentum on faces after advection: "<<sim.Get_Total_Momentum_On_Faces()<<std::endl;
        sim.Identify_Dirichlet();
        T Timing_dirichlet=TIMING_GET;
        sim.Identify_Neumann();
        
        for(int x=3;x<sim.mac_grid.counts.x-3;x++){
            for(int y=2;y<=4;y++){
                sim.cell_neumann(TV_INT(x,y))=true;
                sim.cell_dirichlet(TV_INT(x,y))=false;
                sim.neumann_cell_normal_axis(TV_INT(x,y))=2;}}
        for(int x=3;x<sim.mac_grid.counts.x-3;x++){
            for(int y=sim.mac_grid.counts.y-5;y<=sim.mac_grid.counts.y-3;y++){
                sim.cell_neumann(TV_INT(x,y))=true;
                sim.cell_dirichlet(TV_INT(x,y))=false;
                sim.neumann_cell_normal_axis(TV_INT(x,y))=-2;}}
        for(int x=2;x<=4;x++){
            for(int y=3;y<sim.mac_grid.counts.y-3;y++){
                sim.cell_neumann(TV_INT(x,y))=true;
                sim.cell_dirichlet(TV_INT(x,y))=false;
                sim.neumann_cell_normal_axis(TV_INT(x,y))=1;}}
        for(int x=sim.mac_grid.counts.x-5;x<=sim.mac_grid.counts.x-3;x++){
            for(int y=3;y<sim.mac_grid.counts.y-3;y++){
                sim.cell_neumann(TV_INT(x,y))=true;
                sim.cell_dirichlet(TV_INT(x,y))=false;
                sim.neumann_cell_normal_axis(TV_INT(x,y))=-1;}}
        
        T Timing_neumann=TIMING_GET;
        sim.Classify_Cells();
        T Timing_classify=TIMING_GET;
        sim.Build_Velocity_Divergence();
        T Timing_build_velocity_div=TIMING_GET;
        sim.Solve_For_Pressure();
        T Timing_solve=TIMING_GET;
        sim.Do_Projection();
        T Timing_proj=TIMING_GET;
        // LOG::cout<<"Total momentum on faces after projection: "<<sim.Get_Total_Momentum_On_Faces()<<std::endl;
        sim.Update_Particle_Velocities();
        T Timing_update_v=TIMING_GET;
        // LOG::cout<<"Total momentum on particles after interpolating back: "<<sim.Get_Total_Momentum_On_Particles()<<std::endl;
        // sim.Particle_Based_Body_Collisions();
        sim.Update_Particle_Positions();
        T Timing_update_p=TIMING_GET;
        sim.frame++;
        T Timing_total=Timing_reinit+Timing_weights+Timing_rasterize+Timing_advection+Timing_dirichlet+Timing_neumann+Timing_classify+Timing_build_velocity_div+Timing_solve+Timing_proj+Timing_update_v+Timing_update_p;
        
        LOG::cout<<std::setw(20)<<"reinit "<<std::setw(10)<<Timing_reinit<<std::setw(10)<<Timing_reinit/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"weights "<<std::setw(10)<<Timing_weights<<std::setw(10)<<Timing_weights/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"rasterize"<<std::setw(10)<<Timing_rasterize<<std::setw(10)<<Timing_rasterize/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"advection "<<std::setw(10)<<Timing_advection<<std::setw(10)<<Timing_advection/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"dirichlet "<<std::setw(10)<<Timing_dirichlet<<std::setw(10)<<Timing_dirichlet/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"neumann "<<std::setw(10)<<Timing_neumann<<std::setw(10)<<Timing_neumann/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"classify "<<std::setw(10)<<Timing_classify<<std::setw(10)<<Timing_classify/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"build_velocity_div "<<std::setw(10)<<Timing_build_velocity_div<<std::setw(10)<<Timing_build_velocity_div/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"solve "<<std::setw(10)<<Timing_solve<<std::setw(10)<<Timing_solve/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"proj "<<std::setw(10)<<Timing_proj<<std::setw(10)<<Timing_proj/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"update_v "<<std::setw(10)<<Timing_update_v<<std::setw(10)<<Timing_update_v/Timing_total*100.0<<"%"<<std::endl;
        LOG::cout<<std::setw(20)<<"update_p "<<std::setw(10)<<Timing_update_p<<std::setw(10)<<Timing_update_p/Timing_total*100.0<<"%"<<std::endl;

        total_time_for_solve+=Timing_solve;
        LOG::cout<<"Total time spent in this simulation for CG: "<<total_time_for_solve<<std::endl;

        if(f%frame_jump==0){
            for(int i=0;i<sim.particles.X.m;i++){
                if(i>=air)
                    Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(39.0,64.0,139.0)/255.0);
                else
                    Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,1,1));
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


