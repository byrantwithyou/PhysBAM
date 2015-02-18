#include <Hybrid_Methods/Examples_And_Drivers/MPM_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include "STANDARD_TESTS_2D.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;

    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"run 3D examples");
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;

    parse_args.Parse(true);

    MPM_EXAMPLE<TV>* example=0;
    if(use_3d){}
    else example=new STANDARD_TESTS<VECTOR<T,2> >(stream_type,parse_args);

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    MPM_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
    delete example;

    LOG::Finish_Logging();
    return 0;
}
