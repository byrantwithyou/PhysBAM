//#####################################################################
// Copyright 2008-2009, Jon Gretarsson, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Embedding_Postprocess/EMBEDDING_POSTPROCESS.h"
#include "Muscle_Tests/MUSCLE_TESTS.h"
#include "Particle_Postprocess/PARTICLE_POSTPROCESS.h"
#include "Smoke_Postprocess/SMOKE_POSTPROCESS.h"
#include "Sphere_Postprocess/SPHERE_POSTPROCESS.h"
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Standard_Tests_Water/STANDARD_TESTS_WATER.h"
#include "Swimming_Tests/SWIMMING_TESTS.h"
#include "Swimming_Tests_Water/SWIMMING_TESTS_WATER.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    typedef GRID<TV> T_GRID;
    STREAM_TYPE stream_type((RW()));

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example;

    if(PARSE_ARGS::Find_And_Remove("-swimming_tests",argc,argv)) example=new SWIMMING_TESTS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-swimming_tests_water",argc,argv)) example=new SWIMMING_TESTS_WATER<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-muscle_tests",argc,argv)) example=new MUSCLE_TESTS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-standard_tests_water",argc,argv)) example=new STANDARD_TESTS_WATER<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-smoke_postprocess",argc,argv)) example=new SMOKE_POSTPROCESS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-embedding_postprocess",argc,argv)) example=new EMBEDDING_POSTPROCESS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-particle_postprocess",argc,argv)) example=new PARTICLE_POSTPROCESS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-sphere_postprocess",argc,argv)) example=new SPHERE_POSTPROCESS<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    example->want_mpi_world=true;
    example->Parse(argc,argv);

    if(example->mpi_world->initialized){
        example->solids_fluids_parameters.mpi_solid_fluid=new MPI_SOLID_FLUID<TV>();
        //int xprocs=parse_args->Get_Integer_Value("-xprocs");
        VECTOR<int,3> proc_counts;
        //if(xprocs) proc_counts=VECTOR<int,3>(xprocs,parse_args->Get_Integer_Value("-yprocs"),parse_args->Get_Integer_Value("-zprocs"));
        if(example->solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
            example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<T_GRID>(*example->fluids_parameters.grid,3,false,proc_counts,VECTOR<bool,3>(),
                example->solids_fluids_parameters.mpi_solid_fluid->fluid_group);
            example->solid_body_collection.deformable_body_collection.simulate=false;
            example->solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;}
        else{
            example->fluids_parameters.simulate=false;
            example->solids_fluids_parameters.mpi_solid_fluid->Create_Fluid_Comm_For_Solid_Nodes();}}
    example->Adjust_Output_Directory_For_MPI(example->solids_fluids_parameters.mpi_solid_fluid);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    if(STANDARD_TESTS<T>* standard_tests=dynamic_cast<STANDARD_TESTS<T>*>(example)) standard_tests->Set_Driver(&driver);
    if(STANDARD_TESTS_WATER<T>* standard_tests_water=dynamic_cast<STANDARD_TESTS_WATER<T>*>(example)) standard_tests_water->Set_Driver(&driver);
    driver.Execute_Main_Program();

    delete example;
    return 0;
}
//#####################################################################
