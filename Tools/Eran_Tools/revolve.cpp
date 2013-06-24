#include <Tools/Math_Tools/constants.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <fstream>
#include <iostream>
#include "../../Personal_Libraries/Eran_Library/POLYGONAL_SURFACE.h"

using namespace PhysBAM;
using namespace std;

VECTOR<double,3> rotate(const VECTOR<double,3> &v, double angle, int axis)
{
    int axis1 = (axis)%3+1;
    int axis2 = (axis+1)%3+1;
    double c = cos(angle);
    double s = sin(angle);
    VECTOR<double,3> rotated_v;
    rotated_v[axis] = v[axis];
    rotated_v[axis1] = c*v[axis1] - s*v[axis2];
    rotated_v[axis2] = s*v[axis1] + c*v[axis2];
    return rotated_v;
}

void zero_out_other_components(int axis, VECTOR<double,3> &v)
{
    for (int i = 1; i <= 3; i++) if (axis != i) v[i] = 0;
}

void revolve(const POLYGONAL_SURFACE &input_surface,
             int num_divisions,
             int axis,
             POLYGONAL_SURFACE &revolved_surface,
             bool cap)
{
    int n = input_surface.polygons(1).m;

    int num_vertices = num_divisions*n;
    if (cap) num_vertices += 2;
    revolved_surface.vertices.Resize(num_vertices);

    int num_polygons = n*num_divisions;
    if (cap) num_polygons += num_divisions;
    revolved_surface.polygons.Resize(num_polygons);

    int i;
    double angle = 2*pi / num_divisions;
    for (i = 1; i <= n; i++)
    {
        revolved_surface.vertices(i) = input_surface.vertices(input_surface.polygons(1)(i));
        for (int j = 2; j <= num_divisions; j++)
        {
            revolved_surface.vertices(i + (j-1)*n) = rotate(revolved_surface.vertices(i), 
                                                            (j-1)*angle, axis);
        }
    }

    int cap_start_vtx_idx, cap_end_vtx_idx;
    if (cap)
    {
        cap_start_vtx_idx = num_vertices - 1;
        VECTOR<double,3> start_vtx = revolved_surface.vertices(1);
        zero_out_other_components(axis, start_vtx);
        revolved_surface.vertices(cap_start_vtx_idx) = start_vtx;

        cap_end_vtx_idx = num_vertices;
        VECTOR<double,3> end_vtx = revolved_surface.vertices(n);
        zero_out_other_components(axis, end_vtx);
        revolved_surface.vertices(cap_end_vtx_idx) = end_vtx;
    }

    bool start_on_bottom = (revolved_surface.vertices(1)[axis] < revolved_surface.vertices(n)[axis]);

    for (i = 1; i <= n; i++)
    {
        if (i==n && cap) break;
        for (int j = 1; j <= num_divisions; j++)
        {
            ARRAY<int> &poly = revolved_surface.polygons(j + (i-1)*num_divisions);
            poly.Resize(4);

            int nextj = (j%num_divisions)+1;
            int nexti = (i%n)+1;

            if (cap && start_on_bottom)
            {
                poly(4) = i + (nextj-1)*n;
                poly(3) = i + (j-1)*n;
                poly(2) = nexti + (j-1)*n;
                poly(1) = nexti + (nextj-1)*n;
            }
            else
            {
                poly(1) = i + (nextj-1)*n;
                poly(2) = i + (j-1)*n;
                poly(3) = nexti + (j-1)*n;
                poly(4) = nexti + (nextj-1)*n;
            }
        }
    }

    if (cap)
    {
        int start_idx = (n-1)*num_divisions;

        for (i = 1; i <= num_divisions; i++)
        {
            ARRAY<int> &poly = revolved_surface.polygons(start_idx + i);
            poly.Resize(3);

            int idx1 = 1 + (i-1)*n;
            int idx2 = 1 + (i%num_divisions)*n;
            if (start_on_bottom)
            {
                poly(1) = idx2;
                poly(2) = idx1;
                poly(3) = cap_start_vtx_idx;
            }
            else
            {
                poly(1) = idx1;
                poly(2) = idx2;
                poly(3) = cap_start_vtx_idx;
            }
        }
    }

    if (cap)
    {
        int start_idx = n*num_divisions;

        for (i = 1; i <= num_divisions; i++)
        {
            ARRAY<int> &poly = revolved_surface.polygons(start_idx + i);
            poly.Resize(3);

            int idx1 = n + (i-1)*n;
            int idx2 = n + (i%num_divisions)*n;
            if (start_on_bottom)
            {
                poly(1) = idx1;
                poly(2) = idx2;
                poly(3) = cap_end_vtx_idx;
            }
            else
            {
                poly(1) = idx2;
                poly(2) = idx1;
                poly(3) = cap_end_vtx_idx;
            }
        }
    }
}

int main(int argc, char *argv[])
{
    bool opt_x=false,opt_y=false,opt_z=false,cap=false;
    int num_divisions=4;
    std::string input_filename;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-x",&opt_x, "rotate along x axis (this is the default)");
    parse_args.Add("-y",&opt_y, "rotate along y axis");
    parse_args.Add("-z",&opt_z, "rotate along z axis");
    parse_args.Add("-c",&cap, "cap ends");
    parse_args.Add("-d", &num_divisions, "num div", "number of divisions");
    parse_args.Extra(&input_filename,"ply file", "ply file to revolve");
    parse_args.Parse();
    int axis = 1;

    if (!FILE_UTILITIES::Is_Ply_File(input_filename))
    {
        cerr << "Not a ply file " << input_filename << endl;
        return -1;
    }

    if (opt_x) axis = 1;
    if (opt_y) axis = 2;
    if (opt_z) axis = 3;

    char output_filename[256];
    sprintf(output_filename, "%s_revolve.ply", FILE_UTILITIES::Get_Basename(input_filename).c_str());

    cout << "Input filename: " << input_filename << endl;
    cout << "Output filename: " << output_filename << endl;
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

    POLYGONAL_SURFACE revolved_surface;
    revolve(input_outline, num_divisions, axis, revolved_surface, cap);
    ofstream output(output_filename);
    revolved_surface.Write(output);
}
