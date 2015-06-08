//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Examples_And_Drivers/SOLIDS_DRIVER.h>
//#include "Cutting/CUTTING_EXAMPLE.h"
//#include "Embedded_Circle/EMBEDDED_CIRCLE_EXAMPLE.h"
//#include "Filament/FILAMENT_EXAMPLE.h"
//#include "Mattress/MATTRESS_EXAMPLE.h"
//#include "Melting_Circle/MELTING_CIRCLE.h"
//#include "Rigid_Bodies_Test/RIGID_BODIES_TEST.h"
//#include "Rigid_Particle/RIGID_PARTICLE_EXAMPLE.h"
//#include "Soft_Constraints_Test/SOFT_CONSTRAINTS_TEST.h"
#include "Standard_Tests/STANDARD_TESTS.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;

    PARSE_ARGS parse_args(argc,argv);
    STANDARD_TESTS<T> example(stream_type,parse_args);
    //MATTRESS_EXAMPLE<T> example(stream_type,parse_args);
    //FILAMENT_EXAMPLE<T> example(stream_type,parse_args);
    //EMBEDDED_CIRCLE_EXAMPLE<T> example(stream_type,parse_args);
    //MELTING_CIRCLE<T> example(stream_type,parse_args);
    //SOFT_CONSTRAINTS_TEST<T> example(stream_type,parse_args);
    //RIGID_BODIES_TEST<T> example(stream_type,parse_args);
    //RIGID_PARTICLE_EXAMPLE<T> example(stream_type,parse_args);
    //CUTTING_EXAMPLE<T> example(stream_type,parse_args);
    example.mpi_world=new MPI_WORLD(parse_args);
    example.After_Construction();

    if(example.mpi_world->initialized) example.solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<TV>);
    example.Adjust_Output_Directory_For_MPI(example.solid_body_collection.deformable_body_collection.mpi_solids);

    SOLIDS_DRIVER<TV> driver(example);
    driver.Execute_Main_Program();

    return 0;
}
//#####################################################################
