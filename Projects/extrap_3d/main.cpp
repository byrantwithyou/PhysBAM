//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include <climits>
#include "Standard_Tests/STANDARD_TESTS.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    RW rw=RW();STREAM_TYPE stream_type(rw); // gcc 3.3.2 workaround
    
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example;
    

    example=new STANDARD_TESTS<T>(stream_type);

    example->want_mpi_world=true;
    PARSE_ARGS parse_args(argc,argv);
    example->Parse(parse_args);

    if(example->mpi_world->initialized) example->solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<TV>);
    example->Adjust_Output_Directory_For_MPI(example->solid_body_collection.deformable_body_collection.mpi_solids);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();

    delete example;
    return 0;
}
//#####################################################################
