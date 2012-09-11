#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "WATER_DRIVER.h"
#include "WATER_EXAMPLE.h"

using namespace PhysBAM;

template<class TV> void Execute_Main_Program(STREAM_TYPE& stream_type,PARSE_ARGS& parse_args,MPI_WORLD& mpi_world)
{ 
    typedef VECTOR<int,TV::dimension> TV_INT;

    int threads=1,scale=100;
    parse_args.Add("-scale",&scale,"scale","fine scale grid resolution");
    parse_args.Add("-threads",&threads,"threads","number of threads");
    parse_args.Parse(true);

    WATER_EXAMPLE<TV>* example=new WATER_EXAMPLE<TV>(stream_type,threads);

    RANGE<TV> range(TV(),TV::All_Ones_Vector());TV_INT counts=TV_INT::All_Ones_Vector()*scale;
    example->Initialize_Grid(counts,range);

    example->last_frame=100;
    parse_args.Add("-restart",&example->restart,"frame","restart frame");
    parse_args.Add("-substep",&example->write_substeps_level,"level","output-substep level");
    parse_args.Add("-e",&example->last_frame,"frame","last frame");
    parse_args.Parse();
    example->write_debug_data=true;
    
    TV point1=TV::All_Ones_Vector()*.45,point2=TV::All_Ones_Vector()*.55;point1(2)=0;point2(2)=.05;
    example->source.min_corner=point1;example->source.max_corner=point2;

    if(mpi_world.initialized){
        example->mpi_grid=new MPI_UNIFORM_GRID<GRID<TV> >(example->mac_grid,3);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("/%d",(mpi_world.rank+1));}

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
    
    WATER_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
}

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    PARSE_ARGS parse_args(argc,argv);
    MPI_WORLD mpi_world(parse_args);

    bool opt_3d=false;
    parse_args.Print_Arguments();
    parse_args.Add("-3d",&opt_3d,"run in 3 dimensions");
    parse_args.Parse(true);

    if(opt_3d)
        Execute_Main_Program<VECTOR<T,3> >(stream_type,parse_args,mpi_world);
    else
        Execute_Main_Program<VECTOR<T,2> >(stream_type,parse_args,mpi_world);

    return 0;
}
