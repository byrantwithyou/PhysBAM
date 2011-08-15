//#####################################################################
// Copyright 2007, Avi Robinson-Mosher, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
//#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_DRIVER.h>
//#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "PLS_FSI_DRIVER.h"
#include "PLS_FSI_EXAMPLE.h"
//#include "Standard_Tests/STANDARD_TESTS.h"
#include "SURFACE_TENSION.h"
//using namespace PhysBAM;

int main(int argc,char* argv[])
{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    typedef double T;
#else
    typedef float T;
#endif

    typedef float RW;
    typedef PhysBAM::VECTOR<T,2> TV;
    typedef PhysBAM::GRID<TV> T_GRID;
    PhysBAM::STREAM_TYPE stream_type((RW()));

    PhysBAM::Two_Phase_Flow_2D_Test::PLS_FSI_EXAMPLE<TV>* example=0;
    example=new PhysBAM::Two_Phase_Flow_2D_Test::SURFACE_TENSION<T>(stream_type);
    example->Parse(argc,argv);
    PhysBAM::Two_Phase_Flow_2D_Test::PLS_FSI_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
    delete example;

    return 0;
}
//#####################################################################
