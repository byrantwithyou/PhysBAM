//#####################################################################
// Copyright 2007, Doug Enright, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN  
//##################################################################### 
//
//#####################################################################
// Enright - September 9, 2003
//#####################################################################
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Circle_Example/CIRCLE_EXAMPLE.h"
#include "Incompressible_Drop/INCOMPRESSIBLE_DROP.h"
#include "Oblique_Sod_ST/OBLIQUE_SOD_ST.h"
#include "Sod_ST_2D/SOD_ST_2D.h"
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Wind_Tunnel/WIND_TUNNEL.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
//#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
 //   typedef double T;
  //  typedef float RW;
//#else
    typedef float T;
    typedef float RW;
//#endif
    typedef VECTOR<T,2> TV;

    STREAM_TYPE stream_type((RW()));

    bool incompressible=false;
    if(PARSE_ARGS::Find_And_Remove("-incompressible",argc,argv)) incompressible=true;

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example=0;
    if(PARSE_ARGS::Find_And_Remove("-sod",argc,argv)) example=new SOD_ST_2D<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-oblique",argc,argv)) example=new OBLIQUE_SOD_ST<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-circle",argc,argv)) example=new CIRCLE_EXAMPLE<T>(stream_type,incompressible);
    else if(PARSE_ARGS::Find_And_Remove("-tunnel",argc,argv)) example=new WIND_TUNNEL<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-drop",argc,argv)) example=new INCOMPRESSIBLE_DROP<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type); //default
    example->want_mpi_world=true;
    int xprocs=PARSE_ARGS::Find_And_Remove_Integer("-xprocs",argc,argv),
        yprocs=PARSE_ARGS::Find_And_Remove_Integer("-yprocs",argc,argv);
    example->Parse(argc,argv);

    if(example->mpi_world->initialized){
        example->solids_fluids_parameters.mpi_solid_fluid=new MPI_SOLID_FLUID<TV>();
        if(example->solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
            example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(*example->fluids_parameters.grid,3,false,VECTOR<int,2>(xprocs,yprocs),
                VECTOR<bool,2>(),example->solids_fluids_parameters.mpi_solid_fluid->fluid_group);
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
    return 0;
}
//#####################################################################
