#include <Dynamics/Fluids_Color_Driver/PLS_FC_DRIVER.h>
#include <Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#include "FLUIDS_COLOR_2D.h"
#include "FLUIDS_COLOR_3D.h"

using namespace PhysBAM;

template<class TV>
void Run_Test(PARSE_ARGS& parse_args,STREAM_TYPE stream_type)
{
    FLUIDS_COLOR<TV>* example=new FLUIDS_COLOR<TV>(stream_type,parse_args);

    Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",example->restart);

    PLS_FC_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
    delete example;
}

int main(int argc,char *argv[])
{
    typedef double RW;
    STREAM_TYPE stream_type((RW()));

    bool use_3d=false;
    bool use_float=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"run 3D examples");
    parse_args.Add("-float",&use_float,"run using floats");
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;

    parse_args.Parse(true);
    if(use_float){
        if(use_3d) Run_Test<VECTOR<float,3> >(parse_args,stream_type);
        else Run_Test<VECTOR<float,2> >(parse_args,stream_type);}
    else{
        if(use_3d) Run_Test<VECTOR<double,3> >(parse_args,stream_type);
        else Run_Test<VECTOR<double,2> >(parse_args,stream_type);}

    LOG::Finish_Logging();
    return 0;
}
