//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include "OD_DRIVER.h"
#include "OD_EXAMPLE.h"
#include "STANDARD_TESTS_2D.h"
#include "STANDARD_TESTS_3D.h"

using namespace PhysBAM;

template<class TV>
void Run_Test(PARSE_ARGS& parse_args,STREAM_TYPE stream_type)
{
    OD_EXAMPLE<TV>* example=new STANDARD_TESTS<TV>(stream_type,parse_args);

    Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",example->restart);

    OD_SOLVER<TV>* solver=new OD_SOLVER<TV>(*example);
    OD_DRIVER<TV> driver(*example,*solver);
    driver.Execute_Main_Program();
    delete example;
    delete solver;
}

int main(int argc,char *argv[])
{
    bool type_double=true;
    bool output_double=true;
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);

    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add_Not("-float_output",&type_double,"Output floats");
    parse_args.Add("-double_output",&type_double,"Output doubles");
    parse_args.Add("-3d",&use_3d,"run 3D examples");
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;

    parse_args.Parse(true);
    if(output_double){
        typedef double RW;
        STREAM_TYPE stream_type((RW()));
        if(type_double){
            typedef double T;
            if(use_3d) Run_Test<VECTOR<T,3> >(parse_args,stream_type);
            else Run_Test<VECTOR<T,2> >(parse_args,stream_type);}
        else{
            typedef float T;
            if(use_3d) Run_Test<VECTOR<T,3> >(parse_args,stream_type);
            else Run_Test<VECTOR<T,2> >(parse_args,stream_type);}}
    else{
        typedef float RW;
        STREAM_TYPE stream_type((RW()));
        if(type_double){
            LOG::cout<<"WARNING: Running with floats but outputting doubles?"<<std::endl;
            typedef double T;
            if(use_3d) Run_Test<VECTOR<T,3> >(parse_args,stream_type);
            else Run_Test<VECTOR<T,2> >(parse_args,stream_type);}
        else{
            typedef float T;
            if(use_3d) Run_Test<VECTOR<T,3> >(parse_args,stream_type);
            else Run_Test<VECTOR<T,2> >(parse_args,stream_type);}}

    LOG::Finish_Logging();
    return 0;
}
