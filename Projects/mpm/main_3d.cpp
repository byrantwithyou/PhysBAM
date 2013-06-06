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
#include <omp.h>
using namespace PhysBAM;

template<class TV>
void Run_Simulation(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,3> TV_INT;
    typedef float RW;
    MPM_SIMULATION<TV> sim;
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
    bool use_bridson=false;
    bool use_projection=false;
    sim.xi=(T)0;
    sim.PROFILING=false;
    parse_args.Add("-test",&test_number,"test","test number");
    parse_args.Add("-o",&output_directory,&use_output_directory,"o","output directory");
    parse_args.Add("-dt",&sim.dt,"dt","dt");
    parse_args.Add("-fj",&frame_jump,"fj","frame jump");
    parse_args.Add("-dump_matrix",&sim.dump_matrix,"dump linear system");
    parse_args.Add("-test_system",&sim.test_system,"test linear system");
    parse_args.Add("-bridson",&use_bridson,"use bridson mesh reconstruction on the fly");
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
        case 1:{ // all materials together
            sim.grid.Initialize(TV_INT(1.2*grid_res+1,2.2*grid_res+1,1.2*grid_res+1),RANGE<TV>(TV(-0.6,-0.6,-0.6),TV(0.6,1.6,0.6)));
            
            sim.particles.Add_Randomly_Sampled_Object(RANGE<TV>(TV(-0.1,-0.3,-0.1),TV(0.1,1.3,0.1)),particle_exclude_radius);
            int first_count=sim.particles.number;
            sim.particles.Set_Material_Properties(0,first_count,
                                                  (T)200*density_scale/first_count, // mass per particle
                                                  (0)*ym/((T)2*((T)1+pr)), // mu
                                                  (0)*ym*pr/(((T)1+pr)*((T)1-2*pr))); // lambda
            sim.particles.Set_Plasticity(0,first_count,
                                         false,-1000,1.001, // plasticity_yield
                                         false,-1,1); // plasticity_clamp
            sim.particles.Set_Visco_Plasticity(0,first_count,
                                               false,100, // visco_nu
                                               5000, // visco_tau
                                               0); // visco_kappa
            sim.particles.Set_Initial_State(0,first_count,
//                                            MATRIX<T,TV::m>(2,0,0,0,1,0,0,0,0.5), // Fe
                                            MATRIX<T,TV::m>::Identity_Matrix(),
                                            MATRIX<T,TV::m>::Identity_Matrix(), // Fp
                                            TV()); // initial velocity
//            for(int p=0;p<sim.particles.number;p++)
//                sim.particles.X(p)=MATRIX<T,TV::m>(2,0,0,0,1,0,0,0,0.5)*sim.particles.Xm(p);
//            
            sim.use_gravity=true;
            break;}
                 
        default: PHYSBAM_FATAL_ERROR("Missing test");};
    
    sim.Initialize();
    if(abs(sim.grid.dX(0)-sim.grid.dX(1))>(T)1e-10 || abs(sim.grid.dX(0)-sim.grid.dX(2))>(T)1e-10 || abs(sim.grid.dX(2)-sim.grid.dX(1))>(T)1e-10){
        LOG::cout<<"grid not uniform! dx: "<<sim.grid.dX<<std::endl;
        exit(0);}
    
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,output_directory);
    
    // projection init
    MPM_PROJECTION<TV> projection(sim);
    
    // Zhu and Bridson
    SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<TV> recons_bridson;
    GRID<TV> recons_bridson_grid(TV_INT()+(grid_res*2+1),sim.grid.domain);
    ARRAY<T,TV_INT> phi_bridson;
    if(use_bridson){
        recons_bridson.Initialize(particle_exclude_radius*2);
        recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,1);
        phi_bridson.array+=1;}

    // draw MPM particles
    for(int i=0;i<sim.particles.X.m;i++) Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,1,1));
    
    Flush_Frame<TV>("mpm");
    
    for(int f=1;f<2900977;f++){
        TIMING_START;
        LOG::cout<<"MPM TIMESTEP "<<f<<std::endl;
        sim.Build_Weights_And_Grad_Weights();
        sim.Build_Helper_Structures_For_Constitutive_Model();
        sim.Rasterize_Particle_Data_To_The_Grid();
        if(sim.frame==0) sim.Compute_Particle_Volumes_And_Densities();
        sim.Compute_Grid_Forces();
        if(sim.use_gravity) sim.Apply_Gravity_To_Grid_Forces();
        sim.Update_Velocities_On_Grid();
        // sim.Grid_Based_Body_Collisions();
        sim.Solve_The_Linear_System(); // so far sim.node_V is achieved via MPM
        if(use_projection){
            projection.Reinitialize();
            projection.Identify_Dirichlet_Cells();
            projection.Identify_Neumann_Cells();
            for(int x=3;x<projection.mac_grid.counts.x-3;x++){
                for(int y=3;y<projection.mac_grid.counts.y-3;y++){
                    projection.cell_neumann(TV_INT(x,y,3))=true;
                    projection.cell_dirichlet(TV_INT(x,y,3))=false;
                    projection.neumann_cell_normal_axis(TV_INT(x,y,3))=3;}}
            for(int x=3;x<projection.mac_grid.counts.x-3;x++){
                for(int y=3;y<projection.mac_grid.counts.y-3;y++){
                    projection.cell_neumann(TV_INT(x,y,projection.mac_grid.counts.z-4))=true;
                    projection.cell_dirichlet(TV_INT(x,y,projection.mac_grid.counts.z-4))=false;
                    projection.neumann_cell_normal_axis(TV_INT(x,y,projection.mac_grid.counts.z-4))=-3;}}
            for(int x=3;x<projection.mac_grid.counts.x-3;x++){
                for(int z=3;z<projection.mac_grid.counts.z-3;z++){
                    for(int y=1;y<=3;y++){
                        projection.cell_neumann(TV_INT(x,y,z))=true;
                        projection.cell_dirichlet(TV_INT(x,y,z))=false;
                        projection.neumann_cell_normal_axis(TV_INT(x,y,z))=2;}}}
            for(int x=3;x<projection.mac_grid.counts.x-3;x++){
                for(int z=3;z<projection.mac_grid.counts.z-3;z++){
                    projection.cell_neumann(TV_INT(x,projection.mac_grid.counts.y-4,z))=true;
                    projection.cell_dirichlet(TV_INT(x,projection.mac_grid.counts.y-4,z))=false;
                    projection.neumann_cell_normal_axis(TV_INT(x,projection.mac_grid.counts.y-4,z))=-2;}}
            for(int y=3;y<projection.mac_grid.counts.y-3;y++){
                for(int z=3;z<projection.mac_grid.counts.z-3;z++){
                    projection.cell_neumann(TV_INT(3,y,z))=true;
                    projection.cell_dirichlet(TV_INT(3,y,z))=false;
                    projection.neumann_cell_normal_axis(TV_INT(3,y,z))=1;}}
            for(int y=3;y<projection.mac_grid.counts.y-3;y++){
                for(int z=3;z<projection.mac_grid.counts.z-3;z++){
                    projection.cell_neumann(TV_INT(projection.mac_grid.counts.x-4,y,z))=true;
                    projection.cell_dirichlet(TV_INT(projection.mac_grid.counts.x-4,y,z))=false;
                    projection.neumann_cell_normal_axis(TV_INT(projection.mac_grid.counts.x-4,y,z))=-1;}}
            
            projection.Identify_Nodes_Of_Non_Dirichlet_Cells();
            projection.Velocities_Corners_To_Faces_MPM_Style();
            projection.Build_Velocity_Divergence();
            projection.Solve_For_Pressure();
            projection.Do_Projection();
            projection.Velocities_Faces_To_Corners_MPM_Style(); // this step modifies sim.node_V
        }
        sim.Update_Deformation_Gradient();
        sim.Update_Particle_Velocities();
        // sim.Particle_Based_Body_Collisions();
        sim.Update_Particle_Positions();
        sim.Update_Dirichlet_Box_Positions();
        sim.Update_Colliding_Object_Positions();
        sim.frame++;
        
        TIMING_END("Current MPM time step totally");
        
        if(f%frame_jump==0){
            // draw MPM particles
            for(int i=0;i<sim.particles.X.m;i++){
                Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(1,1,1));
            }

//            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),projection.mac_grid.counts));it.Valid();it.Next()){
//                if(projection.cell_neumann(it.index)){
//                    PHYSBAM_ASSERT(!projection.cell_dirichlet(it.index));
//                    Add_Debug_Particle(projection.mac_grid.Center(it.index),VECTOR<T,3>(0,1,0));}}
//            
            // Zhu and Bridson
            if(use_bridson){
                recons_bridson.Build_Scalar_Field(sim.particles.X,recons_bridson_grid,phi_bridson,1);
                Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(1,0,0));
//                LEVELSET<TV> ls(recons_bridson_grid,phi_bridson,0);
//                ls.Fast_Marching_Method();
//                phi_bridson.array+=particle_exclude_radius*3;
//                ls.Fast_Marching_Method();
//                phi_bridson.array-=particle_exclude_radius*2;
//                ls.Fast_Marching_Method();
//                Dump_Levelset(recons_bridson_grid,phi_bridson,VECTOR<T,3>(0,1,0));
            }

            Flush_Frame<TV>("mpm");}
        LOG::cout<<std::endl;}
}

int main(int argc,char *argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse(true);
    Run_Simulation<VECTOR<double,3> >(parse_args);
    LOG::Finish_Logging();
    return 0;
}


