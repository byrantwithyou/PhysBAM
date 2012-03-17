#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Point_Clouds/PARTICLES.h>
#include <PhysBAM_Tools/Point_Clouds/PARTICLES_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class TV>
void Read_Particles(const STREAM_TYPE& stream_type,const std::string& file)
{
    PARTICLES<TV> particles;
    FILE_UTILITIES::Read_From_File(stream_type,file,particles);
    for(int i=0;i<particles.array_collection->Size();i++)
        particles.Print(std::cout,i);
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false;
    int dim=3;

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-float","data should be in float format");
    parse_args.Add_Option_Argument("-double","data should be in double format");
    parse_args.Add_Integer_Argument("-d",3,"dimension");
    parse_args.Set_Extra_Arguments(1,"<file>","<obj file> obj file to convert");
    parse_args.Parse(argc, argv);
    if(parse_args.Is_Value_Set(("-d"))) dim=parse_args.Get_Integer_Value("-d");

    std::string filename=parse_args.Extra_Arg(1);

    if(parse_args.Get_Option_Value("-float")) type_double=false;
    if(parse_args.Get_Option_Value("-double")) type_double=true;
    STREAM_TYPE stream_type(type_double?STREAM_TYPE(0.0):STREAM_TYPE(0.0f));

    if(!type_double && dim==1) Read_Particles<VECTOR<float,1> >(stream_type,filename);
    if(!type_double && dim==2) Read_Particles<VECTOR<float,2> >(stream_type,filename);
    if(!type_double && dim==3) Read_Particles<VECTOR<float,3> >(stream_type,filename);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    if(type_double && dim==1) Read_Particles<VECTOR<double,1> >(stream_type,filename);
    if(type_double && dim==2) Read_Particles<VECTOR<double,2> >(stream_type,filename);
    if(type_double && dim==3) Read_Particles<VECTOR<double,3> >(stream_type,filename);
#endif
}
