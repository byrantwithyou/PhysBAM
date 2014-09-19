#include "FLUID_STRESS.h"
#include "STRESS_DRIVER.h"
#include "STRESS_EXAMPLE.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;

    PARSE_ARGS parse_args(argc,argv);
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;
    FLUID_STRESS<TV>* example=new FLUID_STRESS<TV>(stream_type,parse_args);

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    STRESS_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
    delete example;

    LOG::Finish_Logging();
    return 0;
}
