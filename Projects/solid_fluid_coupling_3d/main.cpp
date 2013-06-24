//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Joints/JOINT_MESH.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Standard_Tests_Water/STANDARD_TESTS_WATER.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    STREAM_TYPE stream_type((RW()));

    bool opt_water=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-water",&opt_water,"Use water test");
#if 0
    int xprocs=0,zprocs=0;
    parse_args.Add("-xprocs",&xprocs,"procs","Processors in x direction");
    parse_args.Add("-zprocs",&zprocs,"procs","Processors in x direction");
#endif
    parse_args.Parse(true);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example=0;
    if(opt_water) example=new STANDARD_TESTS_WATER<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    example->want_mpi_world=true;
    example->Parse(parse_args);

    std::cout<<"mpi world "<<example->mpi_world->initialized<<std::endl;
#if 0
    if(example->mpi_world->initialized){
        example->solids_fluids_parameters.mpi_solid_fluid=new MPI_SOLID_FLUID<TV>();
        VECTOR<int,3> proc_counts;
        if(xprocs) proc_counts=VECTOR<int,3>(xprocs,1,zprocs);
        if(example->solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
            example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(*example->fluids_parameters.grid,3,false,proc_counts,VECTOR<bool,3>(),
                example->solids_fluids_parameters.mpi_solid_fluid->fluid_group);
            example->solid_body_collection.deformable_body_collection.simulate=false;
            example->solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;}
        else{
            example->fluids_parameters.simulate=false;
            example->solids_fluids_parameters.mpi_solid_fluid->Create_Fluid_Comm_For_Solid_Nodes();}}
#endif
    std::cout<<"TAG "<<example->solids_fluids_parameters.mpi_solid_fluid<<std::endl;
    example->Adjust_Output_Directory_For_MPI(example->solids_fluids_parameters.mpi_solid_fluid);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
//#####################################################################
