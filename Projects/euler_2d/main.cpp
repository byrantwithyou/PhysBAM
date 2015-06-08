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
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Circle_Example/CIRCLE_EXAMPLE.h"
#include "Incompressible_Drop/INCOMPRESSIBLE_DROP.h"
#include "Oblique_Sod_ST/OBLIQUE_SOD_ST.h"
#include "Sod_ST_2D/SOD_ST_2D.h"
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Wind_Tunnel/WIND_TUNNEL.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,2> TV;

    STREAM_TYPE stream_type((RW()));

    bool opt_sod=false,opt_oblique=false,opt_circle=false,opt_tunnel=false,opt_drop=false,incompressible=false;
    int xprocs=0,yprocs=0;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-sod",&opt_sod,"Use sod test");
    parse_args.Add("-oblique",&opt_oblique,"Use oblique test");
    parse_args.Add("-circle",&opt_circle,"Use circle test");
    parse_args.Add("-tunnel",&opt_tunnel,"Use tunnel test");
    parse_args.Add("-drop",&opt_drop,"Use drop test");
    parse_args.Add("-incompressible",&incompressible,"use incompressible");
    parse_args.Add("-xprocs",&xprocs,"procs","Processors in x direction");
    parse_args.Add("-yprocs",&yprocs,"procs","Processors in y direction");
    parse_args.Parse(true);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>* example=0;
    if(opt_sod) example=new SOD_ST_2D<T>(stream_type,parse_args);
    else if(opt_oblique) example=new OBLIQUE_SOD_ST<T>(stream_type,parse_args);
    else if(opt_circle) example=new CIRCLE_EXAMPLE<T>(stream_type,parse_args,incompressible);
    else if(opt_tunnel) example=new WIND_TUNNEL<T>(stream_type,parse_args);
    else if(opt_drop) example=new INCOMPRESSIBLE_DROP<T>(stream_type,parse_args);
    else example=new STANDARD_TESTS<T>(stream_type,parse_args); //default
    example->mpi_world=new MPI_WORLD(parse_args);
    example->After_Construction();

    if(example->mpi_world->initialized){
        example->solids_fluids_parameters.mpi_solid_fluid=new MPI_SOLID_FLUID<TV>();
        if(example->solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
            example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<TV>(*example->fluids_parameters.grid,3,false,VECTOR<int,2>(xprocs,yprocs),
                VECTOR<bool,2>(),example->solids_fluids_parameters.mpi_solid_fluid->fluid_group);
            example->solid_body_collection.deformable_body_collection.simulate=false;
            example->solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;}
        else{
            example->fluids_parameters.simulate=false;
            example->solids_fluids_parameters.mpi_solid_fluid->Create_Fluid_Comm_For_Solid_Nodes();}
    }
    example->Adjust_Output_Directory_For_MPI(example->solids_fluids_parameters.mpi_solid_fluid);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<TV> driver(*example);
    driver.Execute_Main_Program();
    
    delete example;
    return 0;
}
//#####################################################################
