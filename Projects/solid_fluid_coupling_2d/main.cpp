//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Standard_Tests_Water/STANDARD_TESTS_WATER.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    STREAM_TYPE stream_type((RW()));

    bool opt_water=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-water",&opt_water,"Use water test");
    parse_args.Parse(true);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>* example=0;
    if(opt_water) example=new STANDARD_TESTS_WATER<T>(stream_type,parse_args);
    else example=new STANDARD_TESTS<T>(stream_type,parse_args);
    example->mpi_world=new MPI_WORLD(parse_args);
    example->After_Construction();

    if(example->mpi_world->initialized){
        example->solids_fluids_parameters.mpi_solid_fluid=new MPI_SOLID_FLUID<TV>();
        if(example->solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
            example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<TV>(*example->fluids_parameters.grid,3,false,VECTOR<int,2>(),VECTOR<bool,2>(),
                example->solids_fluids_parameters.mpi_solid_fluid->fluid_group);
            example->solid_body_collection.deformable_body_collection.simulate=false;
            example->solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;}
        else{
            example->fluids_parameters.simulate=false;
            example->solids_fluids_parameters.mpi_solid_fluid->Create_Fluid_Comm_For_Solid_Nodes();}}
    example->Adjust_Output_Directory_For_MPI(example->solids_fluids_parameters.mpi_solid_fluid);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<TV> driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
//#####################################################################
