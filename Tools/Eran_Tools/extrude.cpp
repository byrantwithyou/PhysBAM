#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>
#include "../../Personal_Libraries/Eran_Library/POLYGONAL_SURFACE.h"

using namespace PhysBAM;
using namespace std;

void zero_out_other_components(int axis, VECTOR<double,3> &v)
{
    for (int i = 1; i <= 3; i++) if (axis != i) v[i] = 0;
}

int vertex_index(int i, int j, int n) { return i + j * n; }

void extrude(const POLYGONAL_SURFACE &input_surface,
             double extrude_length, 
             int num_divisions,
             int axis,
             POLYGONAL_SURFACE &extruded_surface)
{
    int n = input_surface.polygons(1).m;

    int num_vertices = (num_divisions+1)*n;
    extruded_surface.vertices.Resize(num_vertices);

    int num_polygons = 2 + n*num_divisions;
    extruded_surface.polygons.Resize(num_polygons);

    VECTOR<double,3> displacement;
    displacement[axis] = extrude_length / num_divisions;

    int i;
    for (i = 1; i <= n; i++)
    {
        for (int j = 0; j <= num_divisions; j++)
        {
            extruded_surface.vertices(vertex_index(i,j,n)) = input_surface.vertices(input_surface.polygons(1)(i)) + (double)j * displacement;
        }
    }

    for (i = 1; i <= n; i++)
    {
        int nexti = (i%n)+1;
        for (int j = 1; j <= num_divisions; j++)
        {
            ARRAY<int> &poly = extruded_surface.polygons(j + (i-1)*num_divisions);
            poly.Append(vertex_index(nexti, j-1, n));
            poly.Append(vertex_index(i, j-1, n));
            poly.Append(vertex_index(i, j, n));
            poly.Append(vertex_index(nexti, j, n));
        }
    }

    for (i = 1; i <= n; i++)
        extruded_surface.polygons(num_polygons-1).Append(vertex_index(i, 0, n));
    for (i = n; i >= 1; i--)
        extruded_surface.polygons(num_polygons).Append(vertex_index(i, num_divisions, n));
}

int main(int argc, char *argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-l", 10, "extrude length", "extrude length");
    parse_args.Add_Integer_Argument("-d", 4, "num div", "number of divisions");
    parse_args.Add_Option_Argument("-x", "extrude along x axis");
    parse_args.Add_Option_Argument("-y", "extrude along y axis");
    parse_args.Add_Option_Argument("-z", "extrude along z axis (this is the default)");
    parse_args.Set_Extra_Arguments(1, "<ply file>", "<ply file> ply file to extrude");

    int extraarg = parse_args.Parse(argc, argv);

    char input_filename[256];
    int num_divisions;
    int axis = 3;

    if (extraarg < argc)
        strcpy(input_filename, argv[extraarg]);

    if (!FILE_UTILITIES::Is_Ply_File(input_filename))
    {
        cerr << "Not a ply file " << input_filename << endl;
        return -1;
    }

    num_divisions = parse_args.Get_Integer_Value("-d");
    if (parse_args.Get_Option_Value("-x")) axis = 1;
    if (parse_args.Get_Option_Value("-y")) axis = 2;
    if (parse_args.Get_Option_Value("-z")) axis = 3;

    double extrude_length = parse_args.Get_Double_Value("-l");

    char output_filename[256];
    sprintf(output_filename, "%s_extrude.ply", FILE_UTILITIES::Get_Basename(input_filename).c_str());

    cout << "Input filename: " << input_filename << endl;
    cout << "Output filename: " << output_filename << endl;
    cout << "Length = " << extrude_length << endl;
    cout << "Divisions = " << num_divisions << endl;
    cout << "Axis = " << ((axis==1) ? "x" : (axis==2) ? "y" : "z") << endl;
    
    POLYGONAL_SURFACE input_outline;
    ifstream input(input_filename);
    if (!input)
    {
        cerr << "Input file " << input_filename << " does not exist!" << endl;
        return -1;
    }

    input_outline.Read(input);
    if (input_outline.polygons.m != 1)
    {
        cerr << "Expected a single polygon in input" << endl;
        return -1;
    }

    POLYGONAL_SURFACE extruded_surface;;
    extrude(input_outline, extrude_length, num_divisions, axis, extruded_surface);
    ofstream output(output_filename);
    extruded_surface.Write(output);
}
