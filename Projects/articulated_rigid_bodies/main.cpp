//#####################################################################
// Copyright 2004-2007, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
//#include "ARB/ARB_EXAMPLE_3D.h"
//#include "ARB_Skeleton/ARB_SKELETON_EXAMPLE.h"
#include "Bridge/BRIDGE_EXAMPLE.h"
#include "Chains/CHAINS_EXAMPLE.h"
//#include "ChainSwing/CHAIN_SWING_EXAMPLE.h"
//#include "Clusters/CLUSTER_EXAMPLE.h"
#include "Curl/CURL_EXAMPLE.h"
//#include "Horse/HORSE_EXAMPLE.h"
//#include "Ladder/LADDER_EXAMPLE.h"
#include "Magnets/MAGNETS_EXAMPLE.h"
#include "Mesh/MESH_EXAMPLE.h"
//#include "Net/NET_EXAMPLE.h"
//#include "Simple_Muscle/SIMPLE_MUSCLE_EXAMPLE.h"
//#include "Skeleton/SKELETON_EXAMPLE.h"
//#include "Snake/SNAKE_EXAMPLE.h"
#include "Standard_Tests/STANDARD_TESTS.h"
#include "Tank/TANK_EXAMPLE.h"
//#include "TP/TP_EXAMPLE.h"
//#include "Visible_Human/VISIBLE_HUMAN_EXAMPLE.h"

using namespace PhysBAM;

int main(int argc,char** argv)
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,3> TV;

    EXAMPLE<TV>* example=0;
    if(PARSE_ARGS::Find_And_Remove("-tank",argc,argv)) example=new TANK_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-mesh",argc,argv)) example=new MESH_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-chains",argc,argv)) example=new CHAINS_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-bridge",argc,argv)) example=new BRIDGE_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-magnets",argc,argv)) example=new MAGNETS_EXAMPLE<T>(stream_type);
    else if(PARSE_ARGS::Find_And_Remove("-curl",argc,argv)) example=new CURL_EXAMPLE<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    PARSE_ARGS parse_args(argc,argv);
    example->Parse(parse_args);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* solid_fluid_example=dynamic_cast<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >*>(example);
    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*solid_fluid_example);
    driver.Execute_Main_Program();
    delete example;

    return 0;
}
//#####################################################################
