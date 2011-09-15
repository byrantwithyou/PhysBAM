//#####################################################################
// Copyright 2008, Jon Gretarsson, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Standard_Tests_Water/STANDARD_TESTS_WATER.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    typedef double T;
#else
    typedef float T;
#endif
    typedef float RW;
    typedef VECTOR<T,2> TV;
    STREAM_TYPE stream_type((RW()));

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example;

    if(PARSE_ARGS::Find_And_Remove("-standard_tests_water",argc,argv)) example=new STANDARD_TESTS_WATER<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    example->want_mpi_world=true;
    example->Parse(argc,argv);

    if(example->mpi_world->initialized) example->solid_body_collection.deformable_body_collection.Set_Mpi_Solids(new MPI_SOLIDS<TV>);
    example->Adjust_Output_Directory_For_MPI(example->solid_body_collection.deformable_body_collection.mpi_solids);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    if(STANDARD_TESTS_WATER<T>* standard_tests_water=dynamic_cast<STANDARD_TESTS_WATER<T>*>(example)) standard_tests_water->Set_Driver(&driver);
    if(STANDARD_TESTS<T>* standard_tests=dynamic_cast<STANDARD_TESTS<T>*>(example)) standard_tests->Set_Driver(&driver);
    driver.Execute_Main_Program();

    delete example;
    return 0;
}
//#####################################################################
