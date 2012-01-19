#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include "../../Personal_Libraries/Eran_Library/POLYGONAL_SURFACE.h"
#include <string.h> // for gcc 2.96

using namespace std;
using namespace PhysBAM;

template<class T> void Convert(const char *input_filename, const char *output_filename, PARSE_ARGS &parse_args)
{
    bool zero_based_vertices = parse_args.Get_Option_Value("-zero_based");

    ifstream input(input_filename);
    if (!input)
    {
        cerr << "Input file " << input_filename << " does not exist!" << endl;
        exit(1);
    }

    TRIANGULATED_SURFACE<T> *triangulated_surface=0;
    FILE_UTILITIES::Create_From_File<T>(input_filename,triangulated_surface);
    triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber();

    POLYGONAL_SURFACE polygonal_surface;
    polygonal_surface.zero_based_vertices = zero_based_vertices;

    polygonal_surface.vertices.Resize(triangulated_surface->particles.number);
    for (int i=1; i<=triangulated_surface->particles.number; i++)
        polygonal_surface.vertices(i)=triangulated_surface->particles.X(i);

    polygonal_surface.polygons.Resize(triangulated_surface->triangle_mesh.triangles.m);
    for (int i=1; i<=triangulated_surface->triangle_mesh.triangles.m; i++) {
        polygonal_surface.polygons(i).Exact_Resize(3);
        for(int j=0;j<3;j++) polygonal_surface.polygons(i)(j)=triangulated_surface->triangle_mesh.triangles(j,i);
    }

    ofstream output(output_filename);
    polygonal_surface.Write(output);
}

int main(int argc, char *argv[])
{
    bool type_double = false;

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-zero_based", "zero based vertices");
    parse_args.Set_Extra_Arguments(1, "<tri file>", "<tri file> tri file to convert");

    char input_filename[256];
    int extraarg = parse_args.Parse(argc, argv);
    if (extraarg < argc)
        strcpy(input_filename, argv[extraarg]);
    else
    {
        parse_args.Print_Usage(true);
        return -1;
    }

    if (parse_args.Get_Option_Value("-float")) type_double = false;
    if (parse_args.Get_Option_Value("-double")) type_double = true;

    char output_filename[256];
    sprintf(output_filename, "%s.ply", FILE_UTILITIES::Get_Basename(input_filename).c_str());

    cout << "Input filename: " << input_filename << " [" << (type_double?"double":"float") << "]" << endl;
    cout << "Output filename: " << output_filename << endl;

    if(!type_double) Convert<float>(input_filename, output_filename, parse_args);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert<double>(input_filename, output_filename, parse_args);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
