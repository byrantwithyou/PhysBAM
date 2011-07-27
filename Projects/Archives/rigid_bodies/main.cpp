//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Tamar Shinar, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/RIGIDS_DRIVER.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "MPI_Example/MPI_EXAMPLE.h"
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Standard_Tests/STANDARD_TESTS_RIGIDS_ONLY.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    STREAM_TYPE stream_type((RW()));

    EXAMPLE<TV>* example;

    if(PARSE_ARGS::Find_And_Remove("-rigids_only",argc,argv)) example=new STANDARD_TESTS_RIGIDS_ONLY<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-testmpi",argc,argv)) example=new MPI_EXAMPLE<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    example->want_mpi_world=true;
    example->Parse(argc,argv);

    if(RIGIDS_EXAMPLE<TV>* rigids_example=dynamic_cast<RIGIDS_EXAMPLE<TV>*>(example)){
        if(rigids_example->mpi_world->initialized) rigids_example->mpi_rigids=new MPI_RIGIDS<TV>();
        rigids_example->Adjust_Output_Directory_For_MPI(rigids_example->mpi_rigids);
        RIGIDS_DRIVER<TV> driver(*rigids_example);
        driver.Execute_Main_Program();}
    else{
        SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* solid_fluid_example=dynamic_cast<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >*>(example);
        SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*solid_fluid_example);
        driver.Execute_Main_Program();}
    delete example;

    return 0;
}
//#####################################################################
