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
    typedef double T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;

    STANDARD_TESTS<T> example(stream_type);
    //MATTRESS_EXAMPLE<T> example(stream_type);
    //FILAMENT_EXAMPLE<T> example(stream_type);
    //EMBEDDED_CIRCLE_EXAMPLE<T> example(stream_type);
    //MELTING_CIRCLE<T> example(stream_type);
    //SOFT_CONSTRAINTS_TEST<T> example(stream_type);
    //RIGID_BODIES_TEST<T> example(stream_type);
    //RIGID_PARTICLE_EXAMPLE<T> example(stream_type)_rigid;
    //CUTTING_EXAMPLE<T> example(stream_type)_cutting;
    example.want_mpi_world=true;
    example.Parse(argc,argv);

    if(example.mpi_world->initialized) example.solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<TV>);
    example.Adjust_Output_Directory_For_MPI(example.solid_body_collection.deformable_body_collection.mpi_solids);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(example);
    driver.Execute_Main_Program();

    return 0;
}
//#####################################################################
