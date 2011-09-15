#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include "../../Personal_Libraries/Eran_Library/POLYGONAL_SURFACE.h"
#include <string.h> // for gcc 2.96

using namespace std;
using namespace PhysBAM;

template<class T> VECTOR<T,3> get_centroid(const ARRAY<int> indices, const ARRAY<VECTOR<T,3> > &vectors)
{
    VECTOR<T,3> centroid;
    for (int i = 1; i <= indices.m; i++)
        centroid += vectors(indices(i));
    centroid /= indices.m;
    return centroid;
}

template<class T> TRIANGULATED_SURFACE<T> *
create_triangulated_surface(const POLYGONAL_SURFACE &polygonal_surface, bool flip = false, bool centroid_division = false)
{
    int i;

    // Initialize particles
    SOLIDS_PARTICLES<T,VECTOR<T,3> > *particles = new SOLIDS_PARTICLES<T,VECTOR<T,3> >();
    particles->Increase_Array_Size(polygonal_surface.vertices.m);
    for (i = 1; i <= polygonal_surface.vertices.m; i++)
    {
        int idx = particles->Add_Particle();
        particles->X(idx) = polygonal_surface.vertices(i);
    }
    
    // Initialize triangle mesh
    // Triangulate non-triangular polygons (assumes the polygons
    // are planar and doesn't try to triangulate smartly)
    TRIANGLE_MESH *triangle_mesh = new TRIANGLE_MESH();
    triangle_mesh->number_nodes = polygonal_surface.vertices.m;

    int number_triangles = 0;
    for (i = 1; i <= polygonal_surface.polygons.m; i++)
    {
        if (centroid_division && polygonal_surface.polygons(i).m > 4)
            number_triangles += polygonal_surface.polygons(i).m;
        else
            number_triangles += polygonal_surface.polygons(i).m - 2;
    }

    triangle_mesh->triangles.Resize(3, number_triangles);
    int triidx = 0;
    for (i = 1; i <= polygonal_surface.polygons.m; i++)
    {
        if (centroid_division && polygonal_surface.polygons(i).m > 4)
        {
            VECTOR<T,3> centroid = get_centroid(polygonal_surface.polygons(i), polygonal_surface.vertices);
            int new_particle_idx = particles->Add_Particle();
            particles->X(new_particle_idx) = centroid;
            triangle_mesh->number_nodes++;

            for (int j = 1; j <= polygonal_surface.polygons(i).m; j++)
            {
                triidx++;

                int nextj = (j%polygonal_surface.polygons(i).m)+1;

                if (flip)
                {
                    triangle_mesh->triangles(3, triidx) = new_particle_idx;
                    triangle_mesh->triangles(2, triidx) = polygonal_surface.polygons(i)(j);
                    triangle_mesh->triangles(1, triidx) = polygonal_surface.polygons(i)(nextj);
                }
                else
                {
                    triangle_mesh->triangles(1, triidx) = new_particle_idx;
                    triangle_mesh->triangles(2, triidx) = polygonal_surface.polygons(i)(j);
                    triangle_mesh->triangles(3, triidx) = polygonal_surface.polygons(i)(nextj);
                }
            }
        }
        else
        {
            for (int j = 1; j <= polygonal_surface.polygons(i).m - 2; j++)
            {
                triidx++;

                if (flip)
                {
                    triangle_mesh->triangles(3, triidx) = polygonal_surface.polygons(i)(1);
                    triangle_mesh->triangles(2, triidx) = polygonal_surface.polygons(i)(j+1);
                    triangle_mesh->triangles(1, triidx) = polygonal_surface.polygons(i)(j+2);
                }
                else
                {
                    triangle_mesh->triangles(1, triidx) = polygonal_surface.polygons(i)(1);
                    triangle_mesh->triangles(2, triidx) = polygonal_surface.polygons(i)(j+1);
                    triangle_mesh->triangles(3, triidx) = polygonal_surface.polygons(i)(j+2);
                }
            }
        }
    }

    TRIANGULATED_SURFACE<T> *triangulated_surface = 
        new TRIANGULATED_SURFACE<T>(*triangle_mesh, *particles);

    return triangulated_surface;
}

template<class T> void Convert(const std::string &input_filename, const std::string &output_filename, PARSE_ARGS &parse_args)
{
    bool flip = parse_args.Get_Option_Value("-f");
    bool centroid_division = parse_args.Get_Option_Value("-c");
    bool zero_based_vertices = parse_args.Get_Option_Value("-zero_based");

    ifstream input(input_filename.c_str());
    if (!input)
    {
        cerr << "Input file " << input_filename << " does not exist!" << endl;
        exit(1);
    }

    POLYGONAL_SURFACE polygonal_surface;
    polygonal_surface.zero_based_vertices = zero_based_vertices;
    polygonal_surface.Read(input);

    TRIANGULATED_SURFACE<T> *triangulated_surface = 
        create_triangulated_surface<T>(polygonal_surface, flip, centroid_division);

    ofstream output(output_filename.c_str(), std::ios::binary);
    triangulated_surface->template Write<T>(output);
}

int main(int argc, char *argv[])
{
    bool type_double = false;

    PARSE_ARGS parse_args;
    parse_args.Use_Help_Option(true);
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-f", "flip orientation");
    parse_args.Add_Option_Argument("-c", "use centroid division for polygons with more than 4 vertices");
    parse_args.Add_Option_Argument("-zero_based", "zero based vertices");
    parse_args.Add_String_Argument("-o", "", "output filename");
    parse_args.Set_Extra_Arguments(1, "<ply file>", "<ply file> ply file to convert");

    std::string input_filename;
    int extraarg = parse_args.Parse(argc, argv);
    if (extraarg < argc) input_filename = argv[extraarg];
    else parse_args.Print_Usage(true);

    if (parse_args.Get_Option_Value("-float")) type_double = false;
    if (parse_args.Get_Option_Value("-double")) type_double = true;

    std::string output_filename;
    if (parse_args.Is_Value_Set("-o")) output_filename = parse_args.Get_String_Value("-o");
    else output_filename = FILE_UTILITIES::Get_Basename(input_filename) + ".tri";

    cout << "Input filename: " << input_filename << endl;
    cout << "Output filename: " << output_filename << " [" << (type_double?"double":"float") << "]" << endl;
    
    if(!type_double) Convert<float>(input_filename, output_filename, parse_args);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Convert<double>(input_filename, output_filename, parse_args);
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
