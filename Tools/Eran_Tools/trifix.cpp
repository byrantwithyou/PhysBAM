#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace PhysBAM;

PARSE_ARGS parse_args;

template<class T>
void Do_It()
{
    std::string output_filename=parse_args.Get_String_Value("-o");

    TRIANGULATED_SURFACE<T> *surface=0;FILE_UTILITIES::Create_From_File<T>(parse_args.Extra_Arg(0),surface);
    surface->triangle_mesh.number_nodes=surface->particles.number;
    std::cout << "Writing " << output_filename << std::endl;
    FILE_UTILITIES::Write_To_File<T>(output_filename,*surface);
}

int main(int argc,char *argv[])
{
    parse_args.Add_String_Argument("-o","merged.tri");
    parse_args.Set_Extra_Arguments(-1,"<tri file> <tri file> ...");
    parse_args.Parse(argc,argv);
    if(parse_args.Num_Extra_Args()<1) return 1;
    Do_It<float>();
}
