#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void Convert(const std::string& input_filename,const std::string& output_filename)
{
    TRIANGULATED_SURFACE<T> ts;
    ts.Read_Obj(input_filename);
    Write_To_File<RW>(output_filename,ts);
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false,compute_using_doubles=false;
    std::string input_filename,output_filename;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-compute_using_doubles",&compute_using_doubles,"compute_using_doubles");
    parse_args.Extra(&input_filename,"obj file","obj file to convert");
    parse_args.Extra(&output_filename,"tri file","output tri file name");
    parse_args.Parse();

    if(!File_Extension_Matches_Ignoring_Compression_Suffix(output_filename,"tri",false)){
        std::cerr<<"Not a tri file: "<<output_filename<<std::endl;
        return -1;}

    if(!type_double){
        if(compute_using_doubles){
            std::cout<<"COMPUTING USING DOUBLES!"<<std::endl;
            Convert<double,float>(input_filename,output_filename);
        }else{Convert<float,float>(input_filename,output_filename);}}
    else Convert<double,double>(input_filename,output_filename);
}
