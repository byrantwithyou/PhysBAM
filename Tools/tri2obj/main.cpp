#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void Convert(const std::string& input_filename,const std::string& output_filename)
{
//    typedef VECTOR<T,3> TV;

    TRIANGULATED_SURFACE<T>* triangulated_surface=0;FILE_UTILITIES::Create_From_File<RW>(input_filename,triangulated_surface);
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(output_filename,false);

    std::string header("# simple obj file format:\n"
        "#   # vertex at coordinates (x,y,z)\n"
        "#   v x y z\n"
        "#   # triangle with vertices a,b,c\n"
        "#   f a b c\n"
        "#   # vertices are indexed starting from 1\n"
        "\n");
    (*output)<<header;

    for(int p=0;p<triangulated_surface->particles.Size();p++)
        (*output)<<STRING_UTILITIES::string_sprintf("v %lg %lg %lg\n",triangulated_surface->particles.X(p)[0],triangulated_surface->particles.X(p)[1],triangulated_surface->particles.X(p)[2]);

    for(int e=0;e<triangulated_surface->mesh.elements.m;e++)
        (*output)<<STRING_UTILITIES::string_sprintf("f %d %d %d\n",triangulated_surface->mesh.elements(e)[0]+1,triangulated_surface->mesh.elements(e)[1]+1,triangulated_surface->mesh.elements(e)[2]+1);
    delete output;
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
    parse_args.Extra(&input_filename,"tri file","tri file to convert");
    parse_args.Extra(&output_filename,"obj file","output obj file name");
    parse_args.Parse();


    if(!FILE_UTILITIES::Is_Tri_File(input_filename)){
        std::cerr<<"Not a tri file: "<<input_filename<<std::endl;
        return -1;}

    if(!type_double){
        if(compute_using_doubles){
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
            std::cout<<"COMPUTING USING DOUBLES!"<<std::endl;
            Convert<double,float>(input_filename,output_filename);
#else
            std::cerr<<"Double support not enabled."<<std::endl;exit(1);
#endif
        }else{Convert<float,float>(input_filename,output_filename);}}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert<double,double>(input_filename,output_filename);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
