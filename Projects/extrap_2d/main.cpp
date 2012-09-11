//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;

    STANDARD_TESTS<T> example(stream_type);
    example.want_mpi_world=true;
    PARSE_ARGS parse_args(argc,argv);
    example.Parse(parse_args);

    if(example.mpi_world->initialized) example.solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<TV>);
    example.Adjust_Output_Directory_For_MPI(example.solid_body_collection.deformable_body_collection.mpi_solids);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(example);
    driver.Execute_Main_Program();

    return 0;
}
//#####################################################################
