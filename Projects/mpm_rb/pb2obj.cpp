#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace PhysBAM;

std::string directory;
std::string output_directory;
int threads=1;

template<class TV>
void Run(PARSE_ARGS& parse_args,STREAM_TYPE stream_type)
{
    try{
#pragma omp parallel for schedule(dynamic)
        for(int i=0;;i++){
            MPM_PARTICLES<TV> particle;
            Read_From_File(stream_type,LOG::sprintf("%s/%d/mpm_particles.gz",directory,i),particle);
            std::ostream* output=Safe_Open_Output(LOG::sprintf("%s/%d.obj",output_directory,i+1),false);
            for(int p=0;p<particle.X.m;p++){
                *output<<"v ";
                particle.X(p).Write_Raw(*output);
                *output<<"\n";}
            delete output;}}
    catch(...){}
}

int main(int argc,char *argv[])
{
    bool type_double=true;
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-threads",&threads,"threads","Number of threads");
    parse_args.Extra(&directory,"dir","simulation directory");

#ifdef USE_OPENMP
    omp_set_num_threads(threads);
#pragma omp parallel
#pragma omp single
    {
        if(omp_get_num_threads()!=threads) PHYSBAM_FATAL_ERROR();
        LOG::cout<<"Running on "<<threads<<" threads"<<std::endl;
    }
#else
    PHYSBAM_ASSERT(threads==1);
#endif

    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;

    parse_args.Parse();

    output_directory=directory+"/mpm_particle_obj";
    Create_Directory(output_directory);

    if(type_double) Run<VECTOR<double,3> >(parse_args,STREAM_TYPE(0.0));
    else Run<VECTOR<float,3> >(parse_args,STREAM_TYPE(0.0f));

    LOG::Finish_Logging();
    return 0;
}
