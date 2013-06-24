#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class TV> void Print_Info(std::string &filename)
{
    GRID<TV> grid;
    FILE_UTILITIES::Read_From_File<typename TV::SCALAR>(filename,grid);
    std::cout << "counts = " << grid.counts << std::endl;
    std::cout << "domain = " << grid.domain << std::endl;
    std::cout << "dX = " << grid.dX << std::endl;
    std::cout << "MAC_offset = " << grid.MAC_offset << std::endl;
}

int main(int argc, char *argv[])
{
    bool type_double=false,opt_1d=false,opt_2d=false,opt_3d=false;
    int dimension=3;
    std::string filename;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-1d",&opt_1d, "force 1d mode");
    parse_args.Add("-2d",&opt_2d, "force 2d mode");
    parse_args.Add("-3d",&opt_3d, "force 3d mode");
    parse_args.Extra(&filename,"filename","filename");
    parse_args.Parse();

    // By default detection dimension from file extension
    if (opt_3d) dimension = 3;
    if (opt_2d) dimension = 2;
    if (opt_1d) dimension = 1;

    cout << "Filename: " << filename << " [" << ((type_double)?"double":"float") << ", " << dimension << "D]" << endl;

    if (dimension == 3)
    {
        if (type_double) Print_Info<VECTOR<double,3> >(filename);
        else Print_Info<VECTOR<float,3> >(filename);
    }
    else if (dimension == 2)
    {
        if (type_double) Print_Info<VECTOR<double,2> >(filename);
        else Print_Info<VECTOR<float,2> >(filename);
    }
    else if (dimension == 1)
    {
        if (type_double) Print_Info<VECTOR<double,1> >(filename);
        else Print_Info<VECTOR<float,1> >(filename);
    }
}
