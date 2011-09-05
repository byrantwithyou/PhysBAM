#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN  
//##################################################################### 
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_RLE.h>
#include "../water_rle_2d/Deep_Water/DEEP_WATER.h"
#include "../water_rle_2d/Standard_Tests/STANDARD_TESTS.h"
#include "Boat/BOAT.h"
#include "Falling_Drop/FALLING_DROP.h"
#include "River/RIVER.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef RLE_GRID_3D<T> T_GRID;

    parse_args.Add_Integer_Argument("-resolution",1);
    parse_args.Add_Double_Argument("-depth",0);
    parse_args.Add_Double_Argument("-vertical",0);
    parse_args.Add_Option_Argument("-slope");
    parse_args.Add_Vector_2D_Argument("-processes_per_dimension",VECTOR<double,2>(0,0));

    int resolution=parse_args.Get_Integer_Value("-resolution");
    T optical_depth=(T)parse_args.Get_Double_Value("-depth");
    T vertical_refinement_depth=(T)parse_args.Get_Double_Value("-vertical");
    bool enforce_refinement_slope=parse_args.Get_Option_Value("-slope");
    VECTOR<int,2> processes_per_dimension(parse_args.Get_Vector_2D_Value("-processes_per_dimension"));

    SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>* example;
    if(PARSE_ARGS::Find_And_Remove("-deep",argc,argv)) example=new DEEP_WATER<T_GRID>(stream_type,options.Subexample(),resolution,optical_depth);
    else if(PARSE_ARGS::Find_And_Remove("-river",argc,argv))  example=new RIVER<T>(stream_type,resolution,optical_depth);
    else if(PARSE_ARGS::Find_And_Remove("-boat",argc,argv)) example=new BOAT<T>(stream_type,resolution,optical_depth);
    else example=new STANDARD_TESTS<T_GRID>(stream_type,options.Subexample(),resolution,optical_depth,vertical_refinement_depth,enforce_refinement_slope);
    example->want_mpi_world=true;
    example->Parse(argc,argv);

    if(example->mpi_world->initialized)
        example->mpi_grid=new MPI_RLE_GRID<T_GRID>(*example->fluids_parameters.grid,false,processes_per_dimension,example->fluids_parameters.periodic);
    else if(example->fluids_parameters.periodic!=VECTOR<bool,3>()){LOG::cerr<<"Periodic domains require MPI."<<std::endl;exit(1);}

    example->Adjust_Output_Directory_For_MPI(example->mpi_grid);

    SOLIDS_FLUIDS_DRIVER_RLE<T_GRID> driver(*example);
    driver.Execute_Main_Program();

    delete example;
    return 0;
}
//#####################################################################
#else
int main() {return 0;}
#endif
