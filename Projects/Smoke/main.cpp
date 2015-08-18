#include <Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
#include "SMOKE_DRIVER.h"
#include "SMOKE_EXAMPLE.h"
#include "SMOKE_PARTICLES.h"

using namespace PhysBAM;

template<class TV> void Execute_Main_Program(STREAM_TYPE& stream_type,PARSE_ARGS& parse_args,MPI_WORLD& mpi_world)
{ 
    typedef VECTOR<int,TV::dimension> TV_INT;

    // ARGS
    int threads=1,scale=200; // default scale = 100
    parse_args.Add("-threads",&threads,"threads","number of threads");
    parse_args.Add("-scale",&scale,"scale","fine scale grid resolution");
    parse_args.Parse(true);
    SMOKE_EXAMPLE<TV>* example=new SMOKE_EXAMPLE<TV>(stream_type,threads);
    RANGE<TV> range(TV(),TV::All_Ones_Vector()*0.5);range.max_corner(1)=1;
    TV_INT counts=TV_INT::All_Ones_Vector()*scale/2;counts(1)=scale;
    example->Initialize_Grid(counts,range);
    LOG::cout<<"Grid dX "<<example->mac_grid.dX<<std::endl;
    example->last_frame=1000;// default last_frame = 100
    parse_args.Add("-restart",&example->restart,"frame","restart frame");
    parse_args.Add("-substeps",&example->write_substeps_level,"level","output-substep level");
    parse_args.Add("-e",&example->last_frame,"frame","last frame");
    parse_args.Add("-d_div",&example->debug_divergence,"output the max velocity divergence after projection");
    parse_args.Add("-a",&example->alpha,"alpha","buoyancy parameter alpha");
    parse_args.Add("-o",&example->output_directory,"dir","output directory");
    parse_args.Add("-eapic_order",&example->eapic_order,"frame","last frame");
    parse_args.Add("-use_eapic",&example->use_eapic,"use EAPIC");
    parse_args.Parse();
    
    if(example->eapic_order==1) example->Set_Weights(new PARTICLE_GRID_WEIGHTS_SPLINE<TV,1>(example->mac_grid,/*threads*/1));
    else if(example->eapic_order==2) example->Set_Weights(new PARTICLE_GRID_WEIGHTS_SPLINE<TV,2>(example->mac_grid,/*threads*/1));
    else if(example->eapic_order==3) example->Set_Weights(new PARTICLE_GRID_WEIGHTS_SPLINE<TV,3>(example->mac_grid,/*threads*/1));
    else PHYSBAM_FATAL_ERROR();

    // INITIALIZE PARTICLES, FILL THE GRID.
    // TODO: X0, X, V, C, mass(=1)
    // RANDOM_NUMBERS<typename TV::SCALAR> rand;
    if(example->use_eapic){
        ARRAY<TV> samples;
        TV a;
        a(0)=0.25;a(1)=0.25;samples.Append(a);
        a(0)=0.75;a(1)=0.25;samples.Append(a);
        a(0)=0.25;a(1)=0.75;samples.Append(a);
        a(0)=0.75;a(1)=0.75;samples.Append(a);
        for(CELL_ITERATOR<TV> iterator(example->mac_grid,2);iterator.Valid();iterator.Next()){
            TV X;
            for(int t=0;t<4;t++){
                // rand.Fill_Uniform(X,iterator.Bounding_Box());
                X=iterator.Bounding_Box().min_corner+samples(t)*example->mac_grid.dX.Min();
                int p=example->particles.Add_Element();
                example->particles.X(p)=X;
                example->particles.X0(p)=X;
                example->particles.V(p)=TV();
                example->particles.C(p)=MATRIX<typename TV::SCALAR,TV::m>();
                example->particles.mass(p)=(typename TV::SCALAR)1;}}}
    
    // SOURCE
    TV point1=TV::All_Ones_Vector()*.23,point2=TV::All_Ones_Vector()*.27;point1(1)=0.05;point2(1)=.15;
    example->source.min_corner=point1;
    example->source.max_corner=point2;

    // MPI CRAP
    if(mpi_world.initialized){
        LOG::cout<<"ERROR: MPI initialized? Shouldn't reach here."<<std::endl;
        example->mpi_grid=new MPI_UNIFORM_GRID<TV>(example->mac_grid,5);
        if(example->mpi_grid->Number_Of_Processors()>1) example->output_directory+=LOG::sprintf("/%d",(mpi_world.rank+1));}

    // OUTPUT DIR
    FILE_UTILITIES::Create_Directory(example->output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(example->output_directory+"/common/log.txt",false);
    
    // DRIVER MAIN LOOP
    SMOKE_DRIVER<TV> driver(*example);
    driver.Execute_Main_Program();
    delete example->mpi_grid;
    delete example;
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
    if(opt_3d) Execute_Main_Program<VECTOR<T,3> >(stream_type,parse_args,mpi_world);
    else Execute_Main_Program<VECTOR<T,2> >(stream_type,parse_args,mpi_world);
    LOG::Finish_Logging();
    return 0;
}
