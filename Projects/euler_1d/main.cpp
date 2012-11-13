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
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Bang_Bang_ST/BANG_BANG_ST.h"
#include "Smooth_Flow/SMOOTH_FLOW.h"
#include "Sod_ST/SOD_ST.h"
#include "Sod_ST_Drop/SOD_ST_DROP.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef double RW;
    typedef VECTOR<T,1> TV;

    STREAM_TYPE stream_type((RW()));

    bool opt_sod=false,opt_bangbang=false,opt_smoothflow=false,opt_drop=false;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-sod",&opt_sod,"Use sod test");
    parse_args.Add("-bangbang",&opt_bangbang,"Use bangbang test");
    parse_args.Add("-smoothflow",&opt_smoothflow,"Use smoothflow test");
    parse_args.Add("-drop",&opt_drop,"Use drop test");
    parse_args.Parse(true);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example=0;
    if(opt_sod) example=new SOD_ST<T>(stream_type);
    else if(opt_bangbang) example=new BANG_BANG_ST<T>(stream_type);
    else if(opt_smoothflow) example=new SMOOTH_FLOW<T>(stream_type);
    else if(opt_drop) example=new SOD_ST_DROP<T>(stream_type);
    else example=new SOD_ST<T>(stream_type); //default
    example->want_mpi_world=true;
    example->Parse(parse_args);

    if(example->mpi_world->initialized) example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(*example->fluids_parameters.grid,3);
    example->Adjust_Output_Directory_For_MPI(example->fluids_parameters.mpi_grid);
    
    LOG::cout<<std::setprecision(16)<<std::endl;
    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();
    
    delete example;
    return 0;
}
//#####################################################################

