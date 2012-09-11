#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_REFINEMENT_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_REFINEMENT_EXAMPLE.h>
#include "Smoke_Tests/SMOKE_TESTS.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    LOG::Initialize_Logging(false,false,1<<31,true);

    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Print_Arguments();
    MPI_WORLD mpi_world(parse_args);

    bool opt_3d=true;
    parse_args.Add("-3d",&opt_3d,"do 3d solver");
    parse_args.Parse(true);

    if(!opt_3d){
        typedef VECTOR<T,2> TV;
        typedef VECTOR<int,TV::dimension> TV_INT;

        SMOKE_TESTS<TV>* example=new SMOKE_TESTS<TV>(stream_type,parse_args);

        if(mpi_world.initialized){
            example->coarse_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->coarse_mac_grid,3);
            example->fine_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->fine_mac_grid,3);
            if(example->fine_mpi_grid->Number_Of_Processors()>1){
                example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->fine_mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}}
        FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
        
        INCOMPRESSIBLE_REFINEMENT_DRIVER<TV> driver(*example);
        driver.Execute_Main_Program();}
    else{
        typedef VECTOR<T,3> TV;
        typedef VECTOR<int,TV::dimension> TV_INT;

        SMOKE_TESTS<TV>* example=new SMOKE_TESTS<TV>(stream_type,parse_args);

        if(mpi_world.initialized){
            example->coarse_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->coarse_mac_grid,3);
            example->fine_mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->fine_mac_grid,3);
            if(example->fine_mpi_grid->Number_Of_Processors()>1){
                example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->fine_mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}}
        FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
        
        INCOMPRESSIBLE_REFINEMENT_DRIVER<TV> driver(*example);
        driver.Execute_Main_Program();}

    return 0;
}
