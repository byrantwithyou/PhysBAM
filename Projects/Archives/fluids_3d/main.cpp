//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// MAIN
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
//#include "Armadillo_Clothes/ARMADILLO_CLOTHES.h"
//#include "Bouncing_Drop/BOUNCING_DROP.h"
//#include "Boundary_Test/BOUNDARY_TEST.h"
//#include "Burning_Paper/BURNING_PAPER.h"
#include "DSD_Fire_Ball/DSD_FIRE_BALL_EXAMPLE.h"
#include "DSD_No_Navier_Stokes/DSD_NO_NAVIER_STOKES.h"
//#include "Elastic_Drip/ELASTIC_DRIP.h"
//#include "Falling_Drop/FALLING_DROP.h"
//#include "Filling_Box/FILLING_BOX.h"
//#include "Flat_Surface/FLAT_SURFACE.h"
//#include "Flow_Over_Blob/FLOW_OVER_BLOB.h"
//#include "Flow_Over_Moving_Body/FLOW_OVER_MOVING_BODY.h"
//#include "Flow_Over_Moving_Sphere/FLOW_OVER_MOVING_SPHERE.h"
//#include "Flow_Over_Sphere/FLOW_OVER_SPHERE.h"
#include "Flow_Past_Eftychis/FLOW_PAST_EFTYCHIS.h"
//#include "Flow_Past_Sphere/FLOW_PAST_SPHERE.h"
#include "Glass/GLASS.h"
//#include "Head_ILM/HEAD_ILM.h"
#include "Lighthouse/LIGHTHOUSE.h"
#include "Mass_Conservation/MASS_CONSERVATION.h"
//#include "Melting_Test/MELTING_TEST.h"
//#include "Melting_Test_S3D/MELTING_TEST_S3D.h"
#include "Multiphase_Fire_Examples/MULTIPHASE_FIRE_EXAMPLES.h"
//#include "Plume_Vorticity/PLUME_VORTICITY.h"
//#include "Rising_Smoke/RISING_SMOKE.h"
//#include "Rising_Sphere/RISING_SPHERE.h"
//#include "Rotating_Box/ROTATING_BOX.h"
#include "Sheeting/SHEETING.h"
//#include "Simple_Flame/SIMPLE_FLAME.h"
//#include "Solid_Fluid_Coupling_Test/SOLID_FLUID_COUPLING_TEST.h"
//#include "SPHERE_EMITTER/SPHERE_EMITTER.h"
//#include "Sphere_ILM/SPHERE_ILM.h"
//#include "Sphere_Into_Water/SPHERE_INTO_WATER.h"
//#include "Splash/SPLASH.h"
#include "Spray_And_Foam/SPRAY_AND_FOAM.h"
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Standard_Tests_Multiphase/STANDARD_TESTS_MULTIPHASE.h"
#include "Standard_Tests_Smoke/STANDARD_TESTS_SMOKE.h"
#include "Standard_Tests_SPH/STANDARD_TESTS_SPH.h"
//#include "Surface_Tension/SURFACE_TENSION.h"
//#include "Thin_Shells_Test_Smoke_And_Fire/THIN_SHELLS_TEST_SMOKE_AND_FIRE.h"
//#include "Thin_Shells_Test_Water/THIN_SHELLS_TEST_WATER.h"
#include "Two_Phase/TWO_PHASE.h"
//#include "Underwater_Dome/UNDERWATER_DOME.h"
//#include "Viscous_Letters_Multiphase/VISCOUS_LETTERS_MULTIPHASE.h"
//#include "Vortex_Particles/VORTEX_PARTICLES.h"
//#include "Vortex_Particles_Stream/VORTEX_PARTICLES_STREAM.h"
//#include "Vortex_Test/VORTEX_TEST.h"
//#include "Vortex_Test_Water/VORTEX_TEST_WATER.h"
#include "Waterfall/WATERFALL.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;
    typedef float RW;
    typedef VECTOR<T,3> TV;
    STREAM_TYPE stream_type((RW()));

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* example=0;
    if(PARSE_ARGS::Find_And_Remove("sph",argc,argv)) example=new STANDARD_TESTS_SPH<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("multiphase",argc,argv)) example=new STANDARD_TESTS_MULTIPHASE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("fire",argc,argv)) example=new MULTIPHASE_FIRE_EXAMPLES<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("smoke",argc,argv)) example=new STANDARD_TESTS_SMOKE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("glass",argc,argv)) example=new GLASS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("lighthouse",argc,argv)) example=new LIGHTHOUSE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("sheeting",argc,argv)) example=new SHEETING<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("spray",argc,argv)) example=new SPRAY_AND_FOAM<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("mass_conservation",argc,argv)) example=new MASS_CONSERVATION<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("dsd_fire_ball",argc,argv)) example=new DSD_FIRE_BALL_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("dsd_no_navier",argc,argv)) example=new DSD_NO_NAVIER_STOKES<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("twophase",argc,argv)) example=new TWO_PHASE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("efty",argc,argv)) example=new FLOW_PAST_EFTYCHIS<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("waterfall",argc,argv)) example=new WATERFALL<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    //FILLING_BOX<T> example(stream_type);
    //FLOW_OVER_SPHERE<T> example(stream_type);
    //FLOW_OVER_BLOB<T> example(stream_type);
    //SPHERE_INTO_WATER<T> example(stream_type);
    //HEAD_ILM example(stream_type);
    //BOUNDARY_TEST example(stream_type);
    //SPHERE_ILM example(stream_type);
    //FLOW_OVER_MOVING_SPHERE example(stream_type);
    //RISING_SPHERE example(stream_type);
    //SPHERE_EMITTER example(stream_type);
    //FLOW_OVER_MOVING_BODY<T> example(stream_type);
    //BOUNCING_DROP<T> example(stream_type);
    //ELASTIC_DRIP<T> example(stream_type);
    //SURFACE_TENSION<T> example(example_number);
    //FALLING_DROP<T> example(stream_type);
    //MELTING_TEST<T> example(stream_type);
    //SPLASH<T> example(stream_type);
    //VORTEX_PARTICLES<T> example(stream_type);
    //FLAT_SURFACE<T> example(stream_type);
    //PLUME_VORTICITY<T> example(stream_type);
    //THIN_SHELLS_TEST_SMOKE_AND_FIRE<T> example(stream_type);
    //THIN_SHELLS_TEST_WATER<T> example(stream_type);
    //MELTING_TEST<T> example(stream_type);
    //MELTING_TEST_S3D<T> example(stream_type);
    //RISING_SMOKE<T> example(stream_type);
    //ROTATING_BOX<T> example(stream_type);
    //example=new SIMPLE_FLAME<T>(stream_type);
    //example=new VORTEX_TEST<T>(stream_type);
    //example=new VORTEX_PARTICLES_STREAM<T>(stream_type);
    //FLOW_PAST_SPHERE<T> example(stream_type);
    //BURNING_PAPER<T> example(example_number);
    //ARMADILLO_CLOTHES<T> example(stream_type);
    //SPRAY_AND_FOAM<T> example(stream_type);
    //UNDERWATER_DOME<T> example(stream_type);
    //SOLID_FLUID_COUPLING_TEST<T> example(stream_type);
    example->want_mpi_world=true;
    example->Parse(argc,argv);

    if(example->mpi_world->initialized) example->fluids_parameters.mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(*example->fluids_parameters.grid,3);
    example->Adjust_Output_Directory_For_MPI(example->fluids_parameters.mpi_grid);

    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*example);
    driver.Execute_Main_Program();

    delete example;
    return 0;
}
//#####################################################################
