#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Dynamics/PLS_DRIVER.h>
#include <PhysBAM_Dynamics/PLS_EXAMPLE.h>
#include "Water_Tests/WATER_TESTS.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;
//    typedef VECTOR<T,3> TV;

    MPI_WORLD mpi_world(argc,argv);

    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-restart",0,"restart frame");
    parse_args.Add_Integer_Argument("-scale",128,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");

    parse_args.Parse(argc,argv);
    parse_args.Print_Arguments(argc,argv);
    
    WATER_TESTS<TV>* example=new WATER_TESTS<TV>(stream_type,parse_args);

    if(mpi_world.initialized){
        example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("/%d",(mpi_world.rank+1));}

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    PLS_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    return 0;
}
