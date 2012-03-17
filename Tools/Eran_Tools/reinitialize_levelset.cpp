#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <iostream>
#include <string>
#include "../../Public_Library/Level_Sets/OCTREE_LEVELSET.h"
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>

using namespace PhysBAM;
using namespace FILE_UTILITIES;

std::string filename,output_filename;

template<class T> void Reinitialize_Octree(PARSE_ARGS& parse_args)
{
    if(!parse_args.Is_Value_Set("-grid")){
        std::cout << "Need to specify -grid" << std::endl;
        exit(1);
    }

    std::string grid_filename = parse_args.Get_String_Value("-grid");

    ARRAY<T> phi;OCTREE_GRID<T> grid;
    std::cout << "Reading (" << filename << ", " << grid_filename << ")..." << std::flush;
    FILE_UTILITIES::template Read_From_File<T>(grid_filename,grid);
    FILE_UTILITIES::template Read_From_File<T>(filename,phi);
    OCTREE_LEVELSET<T> levelset(grid,phi);
    int band=parse_args.Get_Integer_Value("-band");
    std::cout << "Reinitializing (band=" << band << ")..." << std::flush;
    levelset.Fast_Marching_Method(0,parse_args.Get_Integer_Value("-band")*grid.Get_Minimum_Cell_Size());
    std::cout << "Writing (" << output_filename << ")..." << std::flush;
    FILE_UTILITIES::template Write_To_File<T>(output_filename,phi);
    std::cout << std::endl;
}

int main(int argc, char *argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-uniform");
    parse_args.Add_Option_Argument("-adaptive");
    parse_args.Add_Option_Argument("-2d");
    parse_args.Add_Option_Argument("-3d");
    parse_args.Add_Integer_Argument("-band",0,"half_band_width","half band width for fast marching method (in grid cells)");
    parse_args.Add_String_Argument("-grid","");
    parse_args.Add_String_Argument("-o","");
    parse_args.Set_Extra_Arguments(1, "<filename>");
    parse_args.Parse(argc, argv);
    if(parse_args.Num_Extra_Args()<1) return -1;
    else filename=parse_args.Extra_Arg(1);

    if(parse_args.Is_Value_Set("-o"))
        output_filename=parse_args.Get_String_Value("-o");
    else
        output_filename=FILE_UTILITIES::Get_Basename(filename)+".reinitialized."+FILE_UTILITIES::Get_File_Extension(filename);

    if(parse_args.Get_Option_Value("-adaptive")){
        if(parse_args.Get_Option_Value("-2d")){ // quadtrees
            assert(false);return -1;
        }
        else{ // octrees
            Reinitialize_Octree<float>(parse_args);
        }
    }
    else{
        assert(false);return -1;
    }
}
