//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Dynamics/Coupled_Driver/PLS_FSI_DRIVER.h>
#include <Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "KANG.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    STREAM_TYPE stream_type((RW()));

    PLS_FSI_EXAMPLE<TV>* example=0;
    example=new KANG<T>(stream_type);
    PARSE_ARGS parse_args(argc,argv);
    example->Parse(parse_args);
    PLS_FSI_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
    delete example;

    return 0;
}
//#####################################################################
