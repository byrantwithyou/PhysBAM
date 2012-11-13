//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_DRIVER.h>
#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "KANG.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    typedef GRID<TV> T_GRID;
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
