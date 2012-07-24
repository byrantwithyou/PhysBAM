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
    parse_args.Print_Arguments();
    FLUIDS_COLOR<TV>* example=new FLUIDS_COLOR<TV>(stream_type,parse_args);

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    PLS_FC_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    LOG::Finish_Logging();
    return 0;
}
