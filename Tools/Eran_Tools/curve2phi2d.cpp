#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <fstream>
#include <iostream>
#include "../../Public_Library/Geometry/LEVELSET_IMPLICIT_CURVE.h"

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

bool Divides(double b, double a)
{
    return (ceil(a/b) == a/b);
}

template<class T> GRID_2D<T> Make_Square_Grid(const GRID_2D<T> &input_grid, int boundary_cells)
{
    BOX_2D<T> input_box = input_grid.Domain();

    // Want to make voxels into cubes
    T voxel_size = input_grid.min_dx_dy;
    int adjusted_m, adjusted_n;

    VECTOR<T,2> box_size = input_box.Size();
    VECTOR<T,2> box_center = input_box.Center();

    adjusted_m = (int)floor(box_size.x/voxel_size) + 2*boundary_cells;
    if (!Divides(voxel_size, box_size.x)) adjusted_m += 2;

    adjusted_n = (int)floor(box_size.y/voxel_size) + 2*boundary_cells;
    if (!Divides(voxel_size, box_size.y)) adjusted_n += 2;

    VECTOR<T,2> new_size = voxel_size * VECTOR<T,2>(adjusted_m-1, adjusted_n-1);

    BOX_2D<T> new_box;
    new_box.Reset_Bounds(box_center - (T)0.5 * new_size);
    new_box.Enlarge_To_Include_Point(box_center + (T)0.5 * new_size);

    return GRID_2D<T>(adjusted_m, adjusted_n, new_box);
}

template<class T> void Convert(const char *input_filename, int boundary_cells, PARSE_ARGS &parse_args)
{
    char output_filename[256] = "";
    double grid_size[2];
    int m, n;
    parse_args.Get_Vector2_Value("-g", grid_size);
    m = (int)grid_size[0];
    n = (int)grid_size[1];

    SEGMENTED_CURVE_2D<T> *segmented_curve=0;
    FILE_UTILITIES::Create_From_File<T>(input_filename,segmented_curve);

    segmented_curve->Update_Bounding_Box();
    BOX_2D<T> box = *segmented_curve->bounding_box;
    if (box.Size().x == 0) { box.xmin -= 1; box.xmax += 1; }
    if (box.Size().y == 0) { box.ymin -= 1; box.ymax += 1; }

    // Make a cube grid using suggested m, n, m and boundary cells
    GRID_2D<T> original_grid(m, n, box);
    GRID_2D<T> grid = Make_Square_Grid(original_grid, boundary_cells);

    if (parse_args.Is_Value_Set("-o"))
    {
        strcpy(output_filename, parse_args.Get_String_Value("-o").c_str());
    }
    else
    {
        char basename[256];
        strcpy(basename, Get_Basename(input_filename).c_str());
        sprintf(output_filename, "%s_%dx%d.phi2d", basename, grid.m, grid.n);
    }

    cout << "Input filename: " << input_filename << endl;
    cout << "Output filename: " << output_filename << endl;
    cout << "------------------------------------" << endl;
    cout << "Suggested grid size: m = " << m << ", n = " << n << endl;
    cout << "Ajusted to make cube voxels and added boundary cells (" << boundary_cells << ")..." << endl;
    cout << "New number of cells: m = " << grid.m << ", n = " << grid.n << endl;
    cout << "Voxel size: " << grid.min_dx_dy << endl;
    cout << "------------------------------------" << endl;
    cout << "Original box: " << box << endl;
    cout << "New box: " << grid.Domain() << endl;
    cout << "------------------------------------" << endl;
    cout << endl;

    ARRAYS<VECTOR<T,2> > phi(1, grid.m, 1, grid.n);
    segmented_curve->Calculate_Signed_Distance_Function(grid, phi, true);

    LEVELSET_IMPLICIT_CURVE<T> levelset_implicit_curve(grid, phi);
    
    ofstream output(output_filename, std::ios::binary);
    levelset_implicit_curve.template Write<T>(output);
}

int main(int argc, char *argv[])
{
    bool type_double = false;
    char input_filename[256];
    double grid_size[2] = {50, 50};
    int boundary_cells = 3;

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Vector2_Argument("-g", grid_size, "grid size", "suggested grid size - actual size may be larger");
    parse_args.Add_Integer_Argument("-b", boundary_cells, "boundary cells", "number of cells outside bounding box");
    parse_args.Add_String_Argument("-o", "", "output filename", "output_filename");
    parse_args.Set_Extra_Arguments(1, "<curve2d file>", "<curve2d file> curve2d file to convert");

    int extraarg = parse_args.Parse();

    if (parse_args.Is_Value_Set("-b"))
        boundary_cells = parse_args.Get_Integer_Value("-b");

    if (extraarg < argc)
        strcpy(input_filename, argv[extraarg]);
    else
    {
        parse_args.Print_Usage(true);
        return -1;
    }

    if (parse_args.Get_Option_Value("-float")) type_double = false;
    if (parse_args.Get_Option_Value("-double")) type_double = true;

    if (!Is_Curve2D_File(input_filename))
    {
        cerr << "Not a curve2d file: " << input_filename << endl;
        return -1;
    }

    if(!type_double) Convert<float>(input_filename, boundary_cells, parse_args);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert<double>(input_filename, boundary_cells, parse_args);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
