#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_DRIVER.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#include "FLUIDS_COLOR.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef double T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;

    PARSE_ARGS parse_args(argc,argv);
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;
    FLUIDS_COLOR<TV>* example=new FLUIDS_COLOR<TV>(stream_type,parse_args);

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    PLS_FC_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
    delete example;

    LOG::Finish_Logging();
    return 0;
}
