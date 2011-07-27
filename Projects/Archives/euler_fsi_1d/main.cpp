//#####################################################################
// Copyright 2003-2007, Doug Enright, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN  
//##################################################################### 
//
//#####################################################################
// Enright - September 7, 2003
//#####################################################################
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Bang_Bang_ST/BANG_BANG_ST.h"
#include "Piston/PISTON.h"
#include "Smooth_Flow/SMOOTH_FLOW.h"
#include "Sod_ST/SOD_ST.h"
#include "Sod_ST_Drop/SOD_ST_DROP.h"
#include "Standard_Tests/STANDARD_TESTS.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    typedef double T;
    typedef double RW;
#else
    typedef float T;
    typedef float RW;
#endif
    typedef VECTOR<T,1> TV;

    STREAM_TYPE stream_type((RW()));
    
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example=0;
    if(PARSE_ARGS::Find_And_Remove("-sod",argc,argv)) example=new SOD_ST<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-piston",argc,argv)) example=new PISTON<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-bangbang",argc,argv)) example=new BANG_BANG_ST<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-smoothflow",argc,argv)) example=new SMOOTH_FLOW<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-drop",argc,argv)) example=new SOD_ST_DROP<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    example->want_mpi_world=true;
    int xprocs=PARSE_ARGS::Find_And_Remove_Integer("-xprocs",argc,argv);
    example->Parse(argc,argv);

    if(example->mpi_world->initialized){
        example->solids_fluids_parameters.mpi_solid_fluid=new MPI_SOLID_FLUID<TV>();
        if(example->solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
            example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(*example->fluids_parameters.grid,3,false,VECTOR<int,1>(xprocs),
                VECTOR<bool,1>(),example->solids_fluids_parameters.mpi_solid_fluid->fluid_group);
            example->solid_body_collection.deformable_body_collection.simulate=false;
            example->solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;}
        else{
            example->fluids_parameters.simulate=false;
            example->solids_fluids_parameters.mpi_solid_fluid->Create_Fluid_Comm_For_Solid_Nodes();}}
    example->Adjust_Output_Directory_For_MPI(example->solids_fluids_parameters.mpi_solid_fluid);
    
    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();
    
    delete example;
    return 0;
}
//#####################################################################

