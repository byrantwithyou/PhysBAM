//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_DRIVER.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_STATE.h>
#include "COMMON_DATA.h"

using namespace PhysBAM;

namespace PhysBAM
{
template<class T> void
Setup_Example(SHALLOW_WATER_STATE<VECTOR<T,1> >& st,PARSE_ARGS& parse_args,COMMON_DATA<VECTOR<T,1> >& cd);

template<class T> void
Setup_Example(SHALLOW_WATER_STATE<VECTOR<T,2> >& st,PARSE_ARGS& parse_args,COMMON_DATA<VECTOR<T,2> >& cd);
}

template<class TV>
void Run_Test(PARSE_ARGS& parse_args)
{
    SHALLOW_WATER_STATE<TV> state(STREAM_TYPE((typename TV::SCALAR)0));
    COMMON_DATA<TV> cd;
    Setup_Example(state,parse_args,cd);
    SHALLOW_WATER_DRIVER<TV> driver(state);
    driver.Execute_Main_Program();
}

int main(int argc,char *argv[])
{
    bool type_double=true;
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);

    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-3d",&use_3d,"run 3D examples");
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;

    parse_args.Parse(true);
    if(type_double){
        if(use_3d) Run_Test<VECTOR<double,2> >(parse_args);
        else Run_Test<VECTOR<double,1> >(parse_args);}
    else{
        if(use_3d) Run_Test<VECTOR<float,2> >(parse_args);
        else Run_Test<VECTOR<float,1> >(parse_args);}

    LOG::Finish_Logging();
    return 0;
}
