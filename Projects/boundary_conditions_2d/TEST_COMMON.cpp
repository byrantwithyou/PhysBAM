#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <cmath>
#include <iomanip>
#include "ACCURACY_INFO.h"
#include "BOUNDARY_CONDITIONS.h"
#include "TEST_COMMON.h"

template<class TV> TEST_COMMON<TV>::
TEST_COMMON()
{
}

template<class TV> TEST_COMMON<TV>::
~TEST_COMMON()
{
    LOG::Finish_Logging();
}

template<class TV> void TEST_COMMON<TV>::
Init_1()
{
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::Instance()->xml=false;
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    LOG::cout<<std::setprecision(2*sizeof(T));

    sim.Init_1(parse_args);

    parse_args.Add_String_Argument("-o","output","Output Directory");
}
template<class TV> void TEST_COMMON<TV>::
Init_2(int argc,char** argv)
{
    parse_args.Print_Arguments(argc,argv);
    parse_args.Parse(argc,argv);
    sim.Init_2(parse_args);

    output_directory=parse_args.Get_String_Value("-o");

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
}
template<class TV> void TEST_COMMON<TV>::
Init_3()
{
    sim.Init_3();
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",sim.obj.grid);
}

template struct TEST_COMMON<VECTOR<double,1> >;
template struct TEST_COMMON<VECTOR<double,2> >;
