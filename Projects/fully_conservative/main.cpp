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
    bool opt_smoke=false,opt_velocity=false;
    parse_args.Add("-smoke",&opt_smoke,"smoke tests");
    parse_args.Add("-velocity",&opt_velocity,"velocity tests");
    parse_args.Parse(true);

    INCOMPRESSIBLE_EXAMPLE<TV>* example;
    if(opt_smoke) example=new SMOKE_TESTS<TV>(stream_type,parse_args);
    else if(opt_velocity) example=new VELOCITY_TESTS<TV>(stream_type,parse_args);
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

    bool opt_2d=false,opt_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Print_Arguments();
    parse_args.Add("-2d",&opt_2d,"do 2d solve");
    parse_args.Add("-3d",&opt_3d,"do 3d solve");
    parse_args.Parse(true);

    if(opt_3d)
        Execute_Main_Program<VECTOR<T,2> >(stream_type,parse_args,mpi_world);
    else if(opt_2d)
        Execute_Main_Program<VECTOR<T,2> >(stream_type,parse_args,mpi_world);
    else
        Execute_Main_Program<VECTOR<T,1> >(stream_type,parse_args,mpi_world);

    return 0;
}
