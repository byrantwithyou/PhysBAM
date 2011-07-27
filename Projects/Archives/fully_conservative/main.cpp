#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include "Advection_Tests/ADVECTION_TESTS.h"
#include "Smoke_Tests/SMOKE_TESTS.h"
#include "Velocity_Tests/VELOCITY_TESTS.h"

using namespace PhysBAM;

template<class TV> void Execute_Main_Program(STREAM_TYPE& stream_type,PARSE_ARGS& parse_args,MPI_WORLD& mpi_world)
{ 
    INCOMPRESSIBLE_EXAMPLE<TV>* example;
    if(parse_args.Is_Value_Set("-smoke")) example=new SMOKE_TESTS<TV>(stream_type,parse_args);
    else if(parse_args.Is_Value_Set("-velocity")) example=new VELOCITY_TESTS<TV>(stream_type,parse_args);
    else example=new ADVECTION_TESTS<TV>(stream_type,parse_args);
        
    if(mpi_world.initialized){
        example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",example->mpi_grid->Number_Of_Processors(),(mpi_world.rank+1));}
    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);

    INCOMPRESSIBLE_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
}

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    STREAM_TYPE stream_type((RW()));

    MPI_WORLD mpi_world(argc,argv);

    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-scale",64,"grid resolution");
    parse_args.Add_Integer_Argument("-restart",0,"restart");
    parse_args.Add_Double_Argument("-cfl",0.9,"cfl");
    parse_args.Add_Integer_Argument("-test_number",1,"test number");
    parse_args.Add_Integer_Argument("-substep",-1,"level of writing sub-steps");
    parse_args.Add_Integer_Argument("-order",1,"order to use");
    parse_args.Add_Option_Argument("-2d","do 2d solve");
    parse_args.Add_Option_Argument("-3d","do 3d solve");
    parse_args.Add_Option_Argument("-smoke","smoke tests");
    parse_args.Add_Option_Argument("-velocity","velocity tests");
    parse_args.Add_Option_Argument("-conservative","use conservative advection");
    parse_args.Add_Option_Argument("-eno","use eno advection");
    parse_args.Add_Option_Argument("-energy","conserve energy");
  
    parse_args.Parse(argc,argv);
    parse_args.Print_Arguments(argc,argv);
    
    if(parse_args.Is_Value_Set("-3d")){
        Execute_Main_Program<VECTOR<T,2> >(stream_type,parse_args,mpi_world);}
    else if(parse_args.Is_Value_Set("-2d")){
        Execute_Main_Program<VECTOR<T,2> >(stream_type,parse_args,mpi_world);}
    else{
        Execute_Main_Program<VECTOR<T,1> >(stream_type,parse_args,mpi_world);}

    return 0;
}
