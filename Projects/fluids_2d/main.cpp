//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Duc Nguyen, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
//#include "Bouncing_Drop/BOUNCING_DROP.h" // drop bouncing on the ground
//#include "Breaking_Wave/BREAKING_WAVE.h"
#include "Density_Targetting/DENSITY_TARGETTING.h"
//#include "Elastic_Drip/ELASTIC_DRIP.h"
//#include "Falling_Drop/FALLING_DROP.h" // drop falling into a pool of liquid
//#include "Flow_Over_Moving_Body/FLOW_OVER_MOVING_BODY.h" // drop falling into a pool of liquid
#include "Flow_Past_Circle/FLOW_PAST_CIRCLE.h"
//#include "Fluid_Control/FLUID_CONTROL.h"
#include "Glass/GLASS.h"
//#include "Glass_Of_Water/GLASS_OF_WATER.h" // glass of water
//#include "Ground_Element/GROUND_ELEMENT.h"
//#include "Implicit_Viscosity/IMPLICIT_VISCOSITY.h" // drop falling into a pool of liquid
#include "Mass_Conservation/MASS_CONSERVATION.h"
//#include "Melting_Test/MELTING_TEST.h"
//#include "Merging_Flame/MERGING_FLAME.h"
#include "Multiphase_Fire_Examples/MULTIPHASE_FIRE_EXAMPLES.h"
//#include "Neumann_Pocket_Test/NEUMANN_POCKET_TEST.h"
//#include "Plane_Jump/PLANE_JUMP.h"
//#include "Plume_Vorticity/PLUME_VORTICITY.h"
#include "Refinement/REFINEMENT.h"
//#include "Rising_Bubble/RISING_BUBBLE.h" // rising bubble
//#include "Simple_Flame/SIMPLE_FLAME.h"
//#include "Solid_Fluid_Coupling_Test/SOLID_FLUID_COUPLING_TEST.h"
#include "Spinning_Bar/SPINNING_BAR.h"
//#include "Splash/SPLASH.h"
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Standard_Tests_Multiphase/STANDARD_TESTS_MULTIPHASE.h"
#include "Standard_Tests_Smoke/STANDARD_TESTS_SMOKE.h"
#include "Standard_Tests_SPH/STANDARD_TESTS_SPH.h"
//#include "Thin_Shells_Test_Smoke_And_Fire/THIN_SHELLS_TEST_SMOKE_AND_FIRE.h"
//#include "Thin_Shells_Test_Water/THIN_SHELLS_TEST_WATER.h"
#include "Two_Phase/TWO_PHASE.h"
//#include "Water_Into_Box/WATER_INTO_BOX.h" // water jet falling into a box
//#include "Water_Over_Circle/WATER_OVER_CIRCLE.h" //water splashing over a static circle
//#include "Wave/WAVE.h" // wave
using namespace PhysBAM;
int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    STREAM_TYPE stream_type((RW()));

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example=0;
    if(PARSE_ARGS::Find_And_Remove("-sph",argc,argv)) example=new STANDARD_TESTS_SPH<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-multiphase",argc,argv)) example=new STANDARD_TESTS_MULTIPHASE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-twophase",argc,argv)) example=new TWO_PHASE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-fire",argc,argv)) example=new MULTIPHASE_FIRE_EXAMPLES<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-density",argc,argv)) example=new DENSITY_TARGETTING<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-glass",argc,argv)) example=new GLASS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-smoke",argc,argv)) example=new STANDARD_TESTS_SMOKE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-spinning_bar",argc,argv)) example=new SPINNING_BAR<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-refinement",argc,argv)) example=new REFINEMENT<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-flow_cylinder",argc,argv)) example=new FLOW_PAST_CIRCLE<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    //BREAKING_WAVE<T> example(stream_type,omega,epsilon,depth);;
    //RISING_BUBBLE<T> example(stream_type);
    //IMPLICIT_VISCOSITY<T> example(stream_type);
    //WATER_INTO_BOX<T> example(stream_type);
    //WATER_OVER_CIRCLE<T> example(stream_type);
    //FLOW_OVER_MOVING_BODY<T> example(stream_type);
    //BOUNCING_DROP<T> example(stream_type);
    //MELTING_TEST<T> example(stream_type);
    //SPLASH<T> example(stream_type);
    //FALLING_DROP<T> example(stream_type);
    //PLUME_VORTICITY<T> example(stream_type);
    //GROUND_ELEMENT<T> example(stream_type);
    //NEUMANN_POCKET_TEST<T> example(stream_type);
    //THIN_SHELLS_TEST_WATER<T> example(stream_type,example_number);
    //THIN_SHELLS_TEST_SMOKE_AND_FIRE<T> example(stream_type,example_number);
    //SIMPLE_FLAME<T> example(stream_type,argc,argv);
    //MERGING_FLAME<T> example(stream_type);
    //PLANE_JUMP<T> example(stream_type);
    //SOLID_FLUID_COUPLING_TEST<T> example(stream_type,example_number);
    //example=new ELASTIC_DRIP<T>(stream_type);
    //example=new FLUID_CONTROL<T>(stream_type);
    example->want_mpi_world=true;
    example->Parse(argc,argv);
    
    FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters=example->fluids_parameters;
    if(example->mpi_world->initialized)
        fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(*fluids_parameters.grid,3,false,VECTOR<int,2>(),fluids_parameters.periodic);
    else if(fluids_parameters.periodic!=VECTOR<bool,2>()){LOG::cerr<<"Periodic domains require MPI."<<std::endl;exit(1);}

    example->Adjust_Output_Directory_For_MPI(fluids_parameters.mpi_grid);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();

    delete example;
    return 0;
}
//#####################################################################
