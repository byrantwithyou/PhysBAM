#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T> void Convert(const std::string& input_filename,const std::string& output_filename)
{
//    typedef VECTOR<T,3> TV;

    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=0;
    Create_From_File(input_filename,tetrahedralized_volume);
    std::ostream* output=Safe_Open_Output_Raw(output_filename,false);

    std::string header("# simple obj file format:\n"
        "#   # vertex at coordinates (x,y,z)\n"
        "#   v x y z\n"
        "#   # tet with vertices a,b,c,d\n"
        "#   f a b c d\n"
        "#   # vertices are indexed starting from 1\n"
        "\n");
    (*output)<<header;

    for(int p=0;p<tetrahedralized_volume->particles.Size();p++)
        LOG::fprintf(*output,"v %lg %lg %lg\n",tetrahedralized_volume->particles.X(p)[0],tetrahedralized_volume->particles.X(p)[1],tetrahedralized_volume->particles.X(p)[2]);

    for(int e=0;e<tetrahedralized_volume->mesh.elements.m;e++)
        LOG::fprintf(*output,"f %d %d %d %d\n",tetrahedralized_volume->mesh.elements(e)[0],tetrahedralized_volume->mesh.elements(e)[1],tetrahedralized_volume->mesh.elements(e)[2],tetrahedralized_volume->mesh.elements(e)[3]);
    delete output;
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false;
    std::string input_filename,output_filename;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Extra(&input_filename,"tet file","tet file to convert");
    parse_args.Extra(&output_filename,"obj file","output obj file name");
    parse_args.Parse();

    if(!File_Extension_Matches_Ignoring_Compression_Suffix(output_filename,"tet",false)){
        std::cerr<<"Not a tet file: "<<input_filename<<std::endl;
        return -1;}

    if(!type_double) Convert<float>(input_filename,output_filename);
    else Convert<double>(input_filename,output_filename);
}
