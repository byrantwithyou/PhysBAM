#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Dynamics/Drivers/PLS_EXAMPLE.h>
#include <PhysBAM_Dynamics/Drivers/PLS_REFINEMENT_DRIVER.h>
#include "Water_Tests/WATER_TESTS.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    LOG::Initialize_Logging(false,false,1<<30,true);
    
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
//    typedef VECTOR<T,2> TV;
    typedef VECTOR<T,3> TV;

    PARSE_ARGS parse_args(argc,argv);
    MPI_WORLD mpi_world(parse_args);

    parse_args.Print_Arguments();
    WATER_TESTS<TV>* example=new WATER_TESTS<TV>(stream_type,parse_args);

    if(mpi_world.initialized){
        example->coarse_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->coarse_mac_grid,3);
        example->fine_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->fine_mac_grid,3);
        if(example->fine_mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->fine_mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}
    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    PLS_REFINEMENT_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();

    delete example;

    return 0;
}
