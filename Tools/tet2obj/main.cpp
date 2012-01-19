#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

template<class T,class RW> void Convert(const std::string& input_filename,const std::string& output_filename)
{
    typedef VECTOR<T,3> TV;

    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=0;FILE_UTILITIES::Create_From_File<RW>(input_filename,tetrahedralized_volume);
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(output_filename,false);

    std::string header("# simple obj file format:\n"
        "#   # vertex at coordinates (x,y,z)\n"
        "#   v x y z\n"
        "#   # tet with vertices a,b,c,d\n"
        "#   f a b c d\n"
        "#   # vertices are indexed starting from 1\n"
        "\n");
    (*output)<<header;

    for(int p=1;p<=tetrahedralized_volume->particles.array_collection->Size();p++)
        (*output)<<STRING_UTILITIES::string_sprintf("v %lg %lg %lg\n",tetrahedralized_volume->particles.X(p)[1],tetrahedralized_volume->particles.X(p)[2],tetrahedralized_volume->particles.X(p)[3]);

    for(int e=0;e<tetrahedralized_volume->mesh.elements.m;e++)
        (*output)<<STRING_UTILITIES::string_sprintf("f %d %d %d %d\n",tetrahedralized_volume->mesh.elements(e)[1],tetrahedralized_volume->mesh.elements(e)[2],tetrahedralized_volume->mesh.elements(e)[3],tetrahedralized_volume->mesh.elements(e)[4]);
    delete output;
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    Initialize_Read_Write_General_Structures();

    bool type_double=false;

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-float","data is in float format");
    parse_args.Add_Option_Argument("-double","data is in double format");
    parse_args.Add_Option_Argument("-compute_using_doubles");
    parse_args.Set_Extra_Arguments(1,"<tet file>","<tet file> tet file to convert");
    parse_args.Set_Extra_Arguments(2,"<obj file>","<obj file> output obj file name");

    parse_args.Parse(argc, argv);

    std::string input_filename=parse_args.Extra_Arg(1);
    std::string output_filename=parse_args.Extra_Arg(2);

    if(parse_args.Get_Option_Value("-float")) type_double=false;
    if(parse_args.Get_Option_Value("-double")) type_double=true;

    if(!FILE_UTILITIES::Is_Tet_File(input_filename)){
        std::cerr<<"Not a tet file: "<<input_filename<<std::endl;
        return -1;}

    if(!type_double){
        if(parse_args.Get_Option_Value("-compute_using_doubles")){
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
