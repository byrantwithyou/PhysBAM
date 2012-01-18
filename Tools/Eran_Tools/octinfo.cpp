#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <stdio.h>
#include <string.h>

using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class T> void Print_Info(const std::string &filename,PARSE_ARGS &parse_args)
{
    OCTREE_GRID<T> octree_grid;
    FILE_UTILITIES::Read_From_File<T>(filename,octree_grid);

    for(int i=1;i<=octree_grid.number_of_nodes;i++)
        std::cout << "Node " << i << ": " << octree_grid.Node_Locations()(i) << std::endl;
}

int main(int argc, char *argv[])
{
    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Set_Extra_Arguments(1, "<filename>");
    int extraarg = parse_args.Parse(argc, argv);
    std::string filename;
    if (extraarg < argc)
        filename = argv[extraarg];
    else
        return -1;

    Print_Info<float>(filename,parse_args);
}
