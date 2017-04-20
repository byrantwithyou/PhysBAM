#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void Convert(const std::string& input_filename,const std::string& output_filename)
{
//    typedef VECTOR<T,3> TV;

    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=0;Create_From_File<RW>(input_filename,tetrahedralized_volume);
    std::ostream* output=Safe_Open_Output(output_filename,false);

    std::string header("# simple obj file format:\n"
        "#   # vertex at coordinates (x,y,z)\n"
        "#   v x y z\n"
        "#   # tet with vertices a,b,c,d\n"
        "#   f a b c d\n"
        "#   # vertices are indexed starting from 1\n"
        "\n");
    (*output)<<header;

    for(int p=0;p<tetrahedralized_volume->particles.Size();p++)
        (*output)<<LOG::sprintf("v %lg %lg %lg\n",tetrahedralized_volume->particles.X(p)[0],tetrahedralized_volume->particles.X(p)[1],tetrahedralized_volume->particles.X(p)[2]);

    for(int e=0;e<tetrahedralized_volume->mesh.elements.m;e++)
        (*output)<<LOG::sprintf("f %d %d %d %d\n",tetrahedralized_volume->mesh.elements(e)[0],tetrahedralized_volume->mesh.elements(e)[1],tetrahedralized_volume->mesh.elements(e)[2],tetrahedralized_volume->mesh.elements(e)[3]);
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
    parse_args.Extra(&input_filename,"tet file","tet file to convert");
    parse_args.Extra(&output_filename,"obj file","output obj file name");
    parse_args.Parse();

    if(!Is_Tet_File(input_filename)){
        std::cerr<<"Not a tet file: "<<input_filename<<std::endl;
        return -1;}

    if(!type_double){
        if(compute_using_doubles){
            std::cout<<"COMPUTING USING DOUBLES!"<<std::endl;
            Convert<double,float>(input_filename,output_filename);
        }else{Convert<float,float>(input_filename,output_filename);}}
    else Convert<double,double>(input_filename,output_filename);
}
