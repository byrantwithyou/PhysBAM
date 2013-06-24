#include <Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_DRIVER.h>
#include <Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include "Smoke_Tests/SMOKE_TESTS.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    PARSE_ARGS parse_args(argc,argv);
    MPI_WORLD mpi_world(parse_args);

    bool opt_3d=false;
    parse_args.Print_Arguments();
    parse_args.Add("-3d",&opt_3d,"do 3d solve");
    parse_args.Parse(true);
    
    if(!opt_3d){
        typedef VECTOR<T,2> TV;
        SMOKE_TESTS<TV>* example=new SMOKE_TESTS<TV>(stream_type,parse_args);
        
        if(mpi_world.initialized){
            example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
            if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}
        FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

        INCOMPRESSIBLE_DRIVER<TV> driver(*example);
        driver.Execute_Main_Program();}
    else{
        typedef VECTOR<T,3> TV;
        SMOKE_TESTS<TV>* example=new SMOKE_TESTS<TV>(stream_type,parse_args);
        
        if(mpi_world.initialized){
            example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
            if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}
        FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

        INCOMPRESSIBLE_DRIVER<TV> driver(*example);
        driver.Execute_Main_Program();}

    return 0;
}
