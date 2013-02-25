//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
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
