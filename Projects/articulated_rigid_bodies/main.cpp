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

    bool opt_tank=false,opt_mesh=false,opt_chains=false,opt_bridge=false,opt_magnets=false,opt_curl=false;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-tank",&opt_tank,"Use tank test");
    parse_args.Add("-mesh",&opt_mesh,"Use mesh test");
    parse_args.Add("-chains",&opt_chains,"Use chains test");
    parse_args.Add("-bridge",&opt_bridge,"Use bridge test");
    parse_args.Add("-magnets",&opt_magnets,"Use magnets test");
    parse_args.Add("-curl",&opt_curl,"Use curl test");
    parse_args.Parse(true);

    EXAMPLE<TV>* example=0;
    if(opt_tank) example=new TANK_EXAMPLE<T>(stream_type);
    else if(opt_mesh) example=new MESH_EXAMPLE<T>(stream_type);
    else if(opt_chains) example=new CHAINS_EXAMPLE<T>(stream_type);
    else if(opt_bridge) example=new BRIDGE_EXAMPLE<T>(stream_type);
    else if(opt_magnets) example=new MAGNETS_EXAMPLE<T>(stream_type);
    else if(opt_curl) example=new CURL_EXAMPLE<T>(stream_type);
    else example=new STANDARD_TESTS<T>(stream_type);
    example->Parse(parse_args);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >* solid_fluid_example=dynamic_cast<SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >*>(example);
    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> > driver(*solid_fluid_example);
    driver.Execute_Main_Program();
    delete example;

    return 0;
}
//#####################################################################
