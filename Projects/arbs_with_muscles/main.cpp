//#####################################################################
// Copyright 2006-2007, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Surface_Muscle/SURFACE_MUSCLE_EXAMPLE.h"
#include "Visible_Human_Muscle/VISIBLE_HUMAN_MUSCLE_EXAMPLE.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    typedef double T;
#else
    typedef float T;
#endif
    typedef float RW;
    typedef VECTOR<T,3> TV;
    STREAM_TYPE stream_type((RW()));

    PARSE_ARGS parse_args(argc,argv);

//    SURFACE_MUSCLE_EXAMPLE<T> example(stream_type);
    VISIBLE_HUMAN_MUSCLE_EXAMPLE<T> example(stream_type);
    example.want_mpi_world=true;
    example.Parse(parse_args);

    if(example.mpi_world->initialized) example.solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<VECTOR<T,3> >);
    example.Adjust_Output_Directory_For_MPI(example.solid_body_collection.deformable_body_collection.mpi_solids);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(example);
    driver.Execute_Main_Program();

    return 0;
}
