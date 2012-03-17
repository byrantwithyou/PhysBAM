#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void Compute_Errors(const std::string& input_base,const int last_frame,const ARRAY<int>& resolutions)
{
    typedef VECTOR<T,2> TV;

    ARRAY<DEFORMABLE_PARTICLES<TV> > particles(resolutions.m);
    ARRAY<ARRAY<T> > errors(resolutions.m);
    for(int i=particles.m;i>=1;i--) errors(i).Resize(last_frame);

    for(int frame=0;frame<last_frame;frame++){
        for(int i=particles.m;i>=1;i--){
            FILE_UTILITIES::Read_From_File<RW>(input_base+STRING_UTILITIES::string_sprintf("%d/%d/deformable_object_particles",resolutions(i),frame),particles(i));
            errors(i)(frame)=0;
            for(int p=0;p<particles(i).Size();p++)
                errors(i)(frame)+=ARRAYS_COMPUTATIONS::Magnitude_Squared(particles(particles.m).X(p)-particles(i).X(p));}}

    for(int i=particles.m;i>=1;i--){
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(input_base+STRING_UTILITIES::string_sprintf("%d/errors.txt",resolutions(i)),false);
        for(int frame=0;frame<last_frame;frame++)
            (*output)<<errors(i)(frame)<<std::endl;
        delete output;}
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false;

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-float","data is in float format");
    parse_args.Add_Option_Argument("-double","data is in double format");
    parse_args.Add_Option_Argument("-compute_using_doubles");
    parse_args.Add_Integer_Argument("-last_frame",30,"last frame");
    parse_args.Set_Extra_Arguments(1,"<base>","<base> simulation directory path except for resolution");

    ARRAY<int> resolutions;
    int resolution=8;
    for(int i=1;i<=5;i++,resolution*=2)
        resolutions.Append(resolution);
    std::cout<<resolutions<<std::endl;
    

    parse_args.Parse(argc, argv);

    std::string input_filename=parse_args.Extra_Arg(1);

    if(parse_args.Get_Option_Value("-float")) type_double=false;
    if(parse_args.Get_Option_Value("-double")) type_double=true;
    int last_frame=parse_args.Get_Integer_Value("-last_frame");


    if(!type_double){
        if(parse_args.Get_Option_Value("-compute_using_doubles")){
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
            std::cout<<"COMPUTING USING DOUBLES!"<<std::endl;
            Compute_Errors<double,float>(input_filename,last_frame,resolutions);
#else
            std::cerr<<"Double support not enabled."<<std::endl;exit(1);
#endif
        }else{Compute_Errors<float,float>(input_filename,last_frame,resolutions);}}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Compute_Errors<double,double>(input_filename,last_frame,resolutions);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
