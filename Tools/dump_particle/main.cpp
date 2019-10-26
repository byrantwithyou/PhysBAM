#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Particles/PARTICLES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

typedef double T;
typedef VECTOR<T,3> TV;

int main(int argc,char *argv[])
{
    std::string input_filename;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Extra(&input_filename,"file","particle data");
    parse_args.Parse();

    PARTICLES<TV> particles;
    Read_From_File(input_filename,particles);

    for(int p=0;p<particles.number;p++)
    {
        std::cout<<"index = "<<p<<"\n";
        for(ATTRIBUTE_INDEX a(0);a<particles.arrays.m;a++)
            particles.arrays(a)->Print(std::cout,p);
        std::cout<<"\n";
    }

    return 0;
}
