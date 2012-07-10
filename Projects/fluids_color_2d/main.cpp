#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_DRIVER.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#include "FLUIDS_COLOR.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;
//    typedef VECTOR<T,3> TV;

    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-restart",0,"restart frame");
    parse_args.Add_Integer_Argument("-resolution",32,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
    parse_args.Add_Integer_Argument("-last_frame",20,"last simulation frame");

    parse_args.Parse(argc,argv);
    parse_args.Print_Arguments(argc,argv);
    
    FLUIDS_COLOR<TV>* example=new FLUIDS_COLOR<TV>(stream_type,parse_args);

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    PLS_FC_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    return 0;
}
