//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN  
//##################################################################### 
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Incompressible_Drop/INCOMPRESSIBLE_DROP.h"
#include "Sphere_Example/SPHERE_EXAMPLE.h"
#include "Standard_Tests/STANDARD_TESTS.h"

using namespace PhysBAM;

template<class T> void main_program(int argc,char* argv[]){
    typedef T RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,3> TV;

    bool incompressible=false;
    if(PARSE_ARGS::Find_And_Remove("-incompressible",argc,argv)) incompressible=true;

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example=0;
    if(PARSE_ARGS::Find_And_Remove("-sphere",argc,argv)) example=new SPHERE_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-drop",argc,argv)) example=new INCOMPRESSIBLE_DROP<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type,incompressible);
    example->want_mpi_world=true;
    int xprocs=PARSE_ARGS::Find_And_Remove_Integer("-xprocs",argc,argv),
        yprocs=PARSE_ARGS::Find_And_Remove_Integer("-yprocs",argc,argv),
        zprocs=PARSE_ARGS::Find_And_Remove_Integer("-zprocs",argc,argv);
    example->Parse(argc,argv);

    if(example->mpi_world->initialized){
        example->solids_fluids_parameters.mpi_solid_fluid=new MPI_SOLID_FLUID<TV>();
        if(example->solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
            example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(*example->fluids_parameters.grid,3,false,VECTOR<int,3>(xprocs,yprocs,zprocs),
                VECTOR<bool,3>(),example->solids_fluids_parameters.mpi_solid_fluid->fluid_group);
            example->solid_body_collection.deformable_body_collection.simulate=false;
            example->solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;}
        else{
            example->fluids_parameters.simulate=false;
            example->solids_fluids_parameters.mpi_solid_fluid->Create_Fluid_Comm_For_Solid_Nodes();}
    }
    example->Adjust_Output_Directory_For_MPI(example->solids_fluids_parameters.mpi_solid_fluid);
    
    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();
    
    delete example;
}

int main(int argc,char* argv[])
{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(PARSE_ARGS::Find_And_Remove("-double",argc,argv))
        main_program<double>(argc,argv);
    else
#endif
        main_program<float>(argc,argv);
    return 0;
}
//#####################################################################

