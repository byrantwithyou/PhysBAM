#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include "Smoke_Tests/SMOKE_TESTS.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    MPI_WORLD mpi_world(argc,argv);

    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-scale",64,"grid resolution");
    parse_args.Add_Integer_Argument("-restart",0,"restart");
    parse_args.Add_Double_Argument("-vc",.06,"vorticity confinement");
    parse_args.Add_Double_Argument("-source_radius",.05,"radius of source");
    parse_args.Add_Double_Argument("-buoyancy",1,"buoyancy constant");
    parse_args.Add_Integer_Argument("-test_number",1,"test number");
    parse_args.Add_Integer_Argument("-substep",-1,"level of writing sub-steps");
    parse_args.Add_Integer_Argument("-upsample",1,"level of refinement");
    parse_args.Add_Option_Argument("-3d","do 3d solve");
    parse_args.Add_Option_Argument("-conservative","use conservative advection");
  
    parse_args.Parse(argc,argv);
    parse_args.Print_Arguments(argc,argv);
    
    if(!parse_args.Is_Value_Set("-3d")){
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
