#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include "STANDARD_TESTS_2D.h"
#include "STANDARD_TESTS_3D.h"

using namespace PhysBAM;

template<class TV>
void Run_Test(PARSE_ARGS& parse_args,STREAM_TYPE stream_type)
{
    MPM_MAC_EXAMPLE<TV>* example=new STANDARD_TESTS<TV>(stream_type,parse_args);

    Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",example->restart);

    MPM_MAC_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
    delete example;
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
        STREAM_TYPE stream_type((double()));
        if(use_3d) Run_Test<VECTOR<double,3> >(parse_args,stream_type);
        else Run_Test<VECTOR<double,2> >(parse_args,stream_type);}
    else{
        STREAM_TYPE stream_type((float()));
        if(use_3d) Run_Test<VECTOR<float,3> >(parse_args,stream_type);
        else Run_Test<VECTOR<float,2> >(parse_args,stream_type);}

    LOG::Finish_Logging();
    return 0;
}
