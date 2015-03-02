#include <Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include "WATER_DRIVER.h"
#include "WATER_EXAMPLE.h"

using namespace PhysBAM;

template<class T>
void Add_Source(WATER_EXAMPLE<VECTOR<T,1> >* example)
{
    PHYSBAM_FATAL_ERROR();
}

template<class T>
void Add_Source(WATER_EXAMPLE<VECTOR<T,2> >* example)
{
    typedef VECTOR<T,2> TV;
    TV point1,point2;RANGE<TV> source;
    //point1=TV::All_Ones_Vector()*(T).5;point1(0)=.4;point1(1)=.95;point2=TV::All_Ones_Vector()*(T).65;point2(0)=.55;point2(1)=1;
    point1=TV::All_Ones_Vector()*(T).5;point1(0)=.95;point1(1)=.6;point2=TV::All_Ones_Vector()*(T).65;point2(0)=1;point2(1)=.75;
    source.min_corner=point1;source.max_corner=point2;
    example->sources.Append(new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(source));
}

template<class T>
void Add_Source(WATER_EXAMPLE<VECTOR<T,3> >* example)
{
    typedef VECTOR<T,3> TV;
    TV point1,point2;CYLINDER<T> source;
    //point1=TV::All_Ones_Vector()*(T).6;point1(0)=.4;point1(1)=.95;point2=TV::All_Ones_Vector()*(T).6;point2(0)=.4;point2(1)=1;
    point1=TV::All_Ones_Vector()*(T).8;point1(0)=.4;point1(2)=.95;point2=TV::All_Ones_Vector()*(T).8;point2(0)=.4;point2(2)=1;
    source.Set_Endpoints(point1,point2);source.radius=.1;
    IMPLICIT_OBJECT<TV>* analytic=new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(source);
    example->sources.Append(analytic);
}

template<class T>
void Add_Sphere(WATER_EXAMPLE<VECTOR<T,1> >* example)
{
    PHYSBAM_FATAL_ERROR();
}

template<class T>
void Add_Sphere(WATER_EXAMPLE<VECTOR<T,2> >* example)
{
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    std::string model_file_name="../../Public_Data/Rigid_Bodies_2D/circle";
    T radius=0.1;
    int rigid_particle_id=example->rigid_body_collection.Add_Rigid_Body(stream_type,model_file_name,radius,true,true,true,true);
    example->rigid_body_collection.rigid_body_particles.rigid_body(rigid_particle_id)->is_static=false;
}

template<class T>
void Add_Sphere(WATER_EXAMPLE<VECTOR<T,3> >* example)
{
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    std::string model_file_name="../../Public_Data/Rigid_Bodies/sphere";
    T radius=0.1;
    int rigid_particle_id=example->rigid_body_collection.Add_Rigid_Body(stream_type,model_file_name,radius,true,true,true,true);
    example->rigid_body_collection.rigid_body_particles.rigid_body(rigid_particle_id)->is_static=false;
}

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;
    //typedef VECTOR<T,3> TV;
    typedef VECTOR<int,TV::dimension> TV_INT;

    PARSE_ARGS parse_args(argc,argv);
    MPI_WORLD mpi_world(parse_args);

    int scale=128,threads=1;
    parse_args.Print_Arguments();
    parse_args.Add("-scale",&scale,"scale","fine scale grid resolution");
    parse_args.Add("-threads",&threads,"threads","number of threads");
    parse_args.Parse(true);

    WATER_EXAMPLE<TV>* example=new WATER_EXAMPLE<TV>(stream_type,threads);

    example->last_frame=100;
    example->cfl=1;
    parse_args.Add("-restart",&example->restart,"frame","restart frame");
    parse_args.Add("-substeps",&example->write_substeps_level,"level","output-substep level");
    parse_args.Add("-e",&example->last_frame,"frame","last frame");
    parse_args.Add("-cfl",&example->cfl,"cfl","cfl number");
    parse_args.Parse();

    example->Initialize_Grid(TV_INT::All_Ones_Vector()*scale,RANGE<TV>(TV(),TV::All_Ones_Vector()));
    Add_Source(example);
    Add_Sphere(example);

    if(mpi_world.initialized){
        example->mpi_grid=new MPI_UNIFORM_GRID<TV>(example->mac_grid,3);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=STRING_UTILITIES::string_sprintf("/%d",(mpi_world.rank+1));}

    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
    
    WATER_DRIVER<TV> driver(*example);

    driver.Execute_Main_Program();

    return 0;
}
