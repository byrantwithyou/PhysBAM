#include <Hybrid_Methods/Examples_And_Drivers/MPM_DRIVER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include "STANDARD_TESTS_2D.h"
#include "STANDARD_TESTS_3D.h"

using namespace PhysBAM;

template<class TV>
void Run_Test(PARSE_ARGS& parse_args,STREAM_TYPE stream_type)
{
    MPM_EXAMPLE<TV>* example=new STANDARD_TESTS<TV>(stream_type,parse_args);

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    MPM_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
    delete example;
}

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    STREAM_TYPE stream_type((RW()));

    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"run 3D examples");
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;

    parse_args.Parse(true);
    if(use_3d) Run_Test<VECTOR<T,3> >(parse_args,stream_type);
    else Run_Test<VECTOR<T,2> >(parse_args,stream_type);

    LOG::Finish_Logging();
    return 0;
}
