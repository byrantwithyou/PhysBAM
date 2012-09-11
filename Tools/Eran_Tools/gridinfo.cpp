#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
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
    bool type_double = false;
    int dimension = 3;
    std::string filename;

    PARSE_ARGS parse_args;

    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add_Option_Argument("-1d", "force 1d mode");
    parse_args.Add_Option_Argument("-2d", "force 2d mode");
    parse_args.Add_Option_Argument("-3d", "force 3d mode");
    parse_args.Set_Extra_Arguments(1, "<filename>");

    parse_args.Parse();

    filename=parse_args.Extra_Arg(0);

    // By default detection dimension from file extension
    if (parse_args.Get_Option_Value("-3d")) dimension = 3;
    if (parse_args.Get_Option_Value("-2d")) dimension = 2;
    if (parse_args.Get_Option_Value("-1d")) dimension = 1;

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
