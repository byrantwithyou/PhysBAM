#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_3D.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <iostream>

using namespace PhysBAM;

typedef float T;

void Keep_Only_Largest_Negative_Region(GRID_3D<T>& grid,ARRAYS<VECTOR<T,3> >& phi)
{
    ARRAYS<VECTOR<int,3> > colors(grid);ARRAYS<VECTOR<bool,3> > null_edge_is_blocked(grid,1);
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int k=0;k<grid.mn;k++) if(phi(i,j,k)>0) colors(i,j,k)=-1; // make outside regions uncolorable
    FLOOD_FILL_3D flood_fill;flood_fill.Optimize_Fill_For_Single_Cell_Regions(true);
    int number_of_colors=flood_fill.Flood_Fill(colors,null_edge_is_blocked,null_edge_is_blocked,null_edge_is_blocked);
    ARRAY<int> region_size(number_of_colors);
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int k=0;k<grid.mn;k++) if(colors(i,j,k)>0) region_size(colors(i,j,k))++;
    int max_region_size=ARRAY<int>::max(region_size);
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int k=0;k<grid.mn;k++) if(colors(i,j,k)>0 && region_size(colors(i,j,k))<max_region_size) phi(i,j,k)*=-1;
}

int main(int argc, char *argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-i","","input file");
    parse_args.Add_String_Argument("-o","","output file");
    parse_args.Set_Extra_Arguments(0, "");
    int extraarg = parse_args.Parse(argc, argv);

    std::string input_filename=parse_args.Get_String_Value("-i");
    std::string output_filename=parse_args.Get_String_Value("-o");

    if(input_filename.empty() || output_filename.empty()) parse_args.Print_Usage(true);

    GRID_3D<T> grid;
    ARRAYS<VECTOR<T,3> > phi;
    LEVELSET_3D<GRID_3D<T> > levelset(grid,phi);
    std::cout << "Reading from '" << input_filename << "'" << std::endl;
    FILE_UTILITIES::Read_From_File<T>(input_filename,levelset);

    Keep_Only_Largest_Negative_Region(grid,phi);

    std::cout << "Fast marching..." << std::endl;
    FAST_LEVELSET_3D<T> fast_levelset(grid,phi);
    fast_levelset.Set_Band_Width(5);
    fast_levelset.Fast_Marching_Method();

    std::cout << "Writing to '" << input_filename << "'" << std::endl;
    FILE_UTILITIES::Write_To_File<T>(output_filename,levelset);
}
