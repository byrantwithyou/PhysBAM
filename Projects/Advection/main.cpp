#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "ADVECTION_DRIVER.h"
#include "ADVECTION_EXAMPLE.h"

using namespace PhysBAM;

template<class TV> void Execute_Main_Program(STREAM_TYPE& stream_type,PARSE_ARGS& parse_args,MPI_WORLD& mpi_world)
{
    typedef VECTOR<int,TV::dimension> TV_INT;
    
    ADVECTION_EXAMPLE<TV>* example=new ADVECTION_EXAMPLE<TV>(stream_type,parse_args.Get_Integer_Value("-threads"));

    int scale=parse_args.Get_Integer_Value("-scale");
    RANGE<TV> range(TV(),TV::All_Ones_Vector()*5);TV_INT counts=TV_INT::All_Ones_Vector()*scale;
    example->Initialize_Grid(counts,range);
    example->restart=parse_args.Get_Integer_Value("-restart");
    example->last_frame=parse_args.Get_Integer_Value("-e");
    example->write_substeps_level=parse_args.Get_Integer_Value("-substep");

    if(mpi_world.initialized){
        example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("/%d",(mpi_world.rank+1));}

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
    
    ADVECTION_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
}

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    MPI_WORLD mpi_world(argc,argv);

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Integer_Argument("-restart",0,"restart frame");
    parse_args.Add_Integer_Argument("-scale",100,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
    parse_args.Add_Integer_Argument("-e",200,"last frame");
    parse_args.Add_Integer_Argument("-threads",1,"number of threads");
    parse_args.Add_Option_Argument("-3d","run in 3 dimensions");
    parse_args.Add_Option_Argument("-2d","run in 2 dimensions");

    parse_args.Parse();
    parse_args.Print_Arguments();
    
    if(parse_args.Is_Value_Set("-3d")){
        Execute_Main_Program<VECTOR<T,3> >(stream_type,parse_args,mpi_world);}
    else if(parse_args.Is_Value_Set("-2d")){
        Execute_Main_Program<VECTOR<T,2> >(stream_type,parse_args,mpi_world);}
    else{
        Execute_Main_Program<VECTOR<T,1> >(stream_type,parse_args,mpi_world);}

    return 0;
}
