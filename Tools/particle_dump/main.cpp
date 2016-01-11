#include <Tools/Log/LOG.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Particles/PARTICLES.h>
#include <Tools/Particles/PARTICLES_FORWARD.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class TV>
void Read_Particles(const STREAM_TYPE& stream_type,const std::string& file)
{
    PARTICLES<TV> particles;
    FILE_UTILITIES::Read_From_File(stream_type,file,particles);
    for(int i=0;i<particles.Size();i++)
        particles.Print(std::cout,i);
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false;
    int dim=3;

    std::string filename;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-d",&dim,"value","dimension");
    parse_args.Extra(&filename,"obj file","obj file to convert");
    parse_args.Parse();

    STREAM_TYPE stream_type(type_double?STREAM_TYPE(0.0):STREAM_TYPE(0.0f));

    if(!type_double && dim==1) Read_Particles<VECTOR<float,1> >(stream_type,filename);
    if(!type_double && dim==2) Read_Particles<VECTOR<float,2> >(stream_type,filename);
    if(!type_double && dim==3) Read_Particles<VECTOR<float,3> >(stream_type,filename);
    if(type_double && dim==1) Read_Particles<VECTOR<double,1> >(stream_type,filename);
    if(type_double && dim==2) Read_Particles<VECTOR<double,2> >(stream_type,filename);
    if(type_double && dim==3) Read_Particles<VECTOR<double,3> >(stream_type,filename);
}
