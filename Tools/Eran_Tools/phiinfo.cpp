#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Tools/Images/RGB_FILE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <fstream>
#include <iostream>
#include "../../Public_Library/Geometry/LEVELSET_IMPLICIT_SURFACE.h"
#include <stdio.h>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class T> T Max_Magnitude_Value(const ARRAYS<VECTOR<T,3> > &phi)
{
    T max = 0;
    for (int i = phi.m_start; i <= phi.m_end; i++)
        for (int j = phi.n_start; j <= phi.n_end; j++)
            for (int k = phi.mn_start; k <= phi.mn_end; k++)
                if (fabs(phi(i, j, k)) > max)
                    max = fabs(phi(i, j, k));
    return max;
}

// direction: 1 = x, 2 = y, 3 = z
template<class T> void Output_Slices(const char *basename, const ARRAYS<VECTOR<T,3> > &phi, 
                                     bool binary_image, int slice_direction)
{
    int slices, slices_start, width, width_start, height, height_start;
    switch (slice_direction)
    {
        case 1:     slices = phi.m; slices_start = phi.m_start;
                    width = phi.mn; width_start = phi.mn_start;
                    height = phi.n; height_start = phi.n_start;
                    break;

        case 2:     slices = phi.n; slices_start = phi.n_start;
                    width = phi.m; width_start = phi.m_start;
                    height = phi.mn; height_start = phi.mn_start;
                    break;

        default:    slices = phi.mn; slices_start = phi.mn_start;
                    width = phi.m; width_start = phi.m_start;
                    height = phi.n; height_start = phi.n_start;
                    break;
    }

    ARRAYS<VECTOR<VECTOR<T,3> ,2> > input_colors(1, width, 1, height);

    T max_value = Max_Magnitude_Value(phi);

    cout << "Writing slices (direction " << slice_direction << ")..." << flush;

    for (int slice = slices_start; slice < slices_start + slices; slice++)
    {
        cout << slice << " " << flush;

        char output_filename[256];
        sprintf(output_filename, "%s_%03d.rgb", basename, slice);

        for (int x = width_start; x < width_start + width; x++)
            for (int y = height_start; y < height_start + height; y++)
            {
                T value;
                switch (slice_direction)
                {
                    case 1:     value = phi(slice, y, x); break;
                    case 2:     value = phi(x, slice, y); break;
                    default:    value = phi(x, y, slice); break;
                }

                VECTOR<T,3> &color = input_colors(x - width_start + 1, y - height_start + 1);

                if (binary_image)
                {
                    color = (value > 0) ? VECTOR<T,3>(0,0,0) : VECTOR<T,3>(1,1,1);
                }
                else
                {
                    value /= max_value; // Normalize value
                    if (value > 0)
                        color = VECTOR<T,3>(0, 0, 1 - value);
                    else if (value < 0)
                        color = VECTOR<T,3>(1 + value, 0, 0);
                    else
                        color = VECTOR<T,3>(0, 1, 0);
                }
            }

        IMAGE<T>::Write(output_filename,input_colors);
    }

    cout << endl;
}

template<class T> void Print_Info_3D(const char *filename,
                                     bool output_slices,
                                     bool binary_image,
                                     int slice_direction,
                                     bool verbose)
{
    LEVELSET_IMPLICIT_SURFACE<T> *implicit_surface=0;
    FILE_UTILITIES::Create_From_File<T>(filename,implicit_surface);

    cout << "implicit_surface.minimum_cell_size = " << implicit_surface->minimum_cell_size << endl;
    GRID_3D<T> &grid = implicit_surface->levelset.grid;
    cout << "implicit_surface.levelset.grid" << endl;
    cout << "grid m = " << grid.m << ", ";
    cout << "n = " << grid.n << ", ";
    cout << "mn = " << grid.mn << endl;
    cout << "(xmin,xmax) = (" << grid.xmin << ", " << grid.xmax << ")" << endl;
    cout << "(ymin,ymax) = (" << grid.ymin << ", " << grid.ymax << ")" << endl;
    cout << "(zmin,zmax) = (" << grid.zmin << ", " << grid.zmax << ")" << endl;
    cout << "dx = " << grid.dx << ", dy = " << grid.dy << ", dz = " << grid.dz << endl;
    cout << "MAC_offset = " << grid.MAC_offset << endl;
    cout << endl;

    ARRAYS<VECTOR<T,3> > &phi = implicit_surface->levelset.phi;
    cout << "implicit_surface.levelset.phi" << endl;
    cout << "phi m = " << phi.m << ", n = " << phi.n << ", mn = " << phi.mn << endl;
    cout << "(m_start,m_end) = (" << phi.m_start << "," << phi.m_end << ")" << endl;
    cout << "(n_start,n_end) = (" << phi.n_start << "," << phi.n_end << ")" << endl;
    cout << "(mn_start,mn_end) = (" << phi.mn_start << "," << phi.mn_end << ")" << endl;

    if (output_slices)
        Output_Slices("slice", phi, binary_image, slice_direction);

    if (verbose)
    {
        for (int k = 1; k <= phi.mn; k++)
        {
            printf("[slice k=%d]\n", k);
            for (int j = phi.n; j >= 1; j--)
            {
                for (int i = 1; i <= phi.m; i++) printf("%8.3f ", phi(i,j,k));
                printf("\n");
            }
            printf("\n");
        }
    }
}

template<class T> void Print_Info_2D(const char *filename, bool verbose)
{
    GRID_2D<T> grid;
    ARRAYS<VECTOR<T,2> > phi;
    LEVELSET_2D<T> levelset(grid, phi);
    FILE_UTILITIES::Read_From_File<T>(filename,levelset);

    cout << "levelset.grid" << endl;
    cout << "grid m = " << grid.m << ", n = " << grid.n << endl;
    cout << "(xmin,xmax) = (" << grid.xmin << ", " << grid.xmax << ")" << endl;
    cout << "(ymin,ymax) = (" << grid.ymin << ", " << grid.ymax << ")" << endl;
    cout << "dx = " << grid.dx << ", dy = " << grid.dy << endl;
    cout << "MAC_offset = " << grid.MAC_offset << endl;
    cout << endl;

    cout << "levelset.phi" << endl;
    cout << "phi m = " << phi.m << ", n = " << phi.n << endl;
    cout << "(m_start,m_end) = (" << phi.m_start << "," << phi.m_end << ")" << endl;
    cout << "(n_start,n_end) = (" << phi.n_start << "," << phi.n_end << ")" << endl;

    if (verbose)
    {
        for (int j = phi.n_end; j >= phi.n_start; j--) 
        {
            for (int i = phi.m_start; i <= phi.m_end; i++) printf("%8.3f ", phi(i,j));
            printf("\n");
        }
    }
}

int main(int argc, char *argv[])
{
    bool type_double = false;
    int dimension = -1;
    char filename[256];

    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-2d", "force 2d mode");
    parse_args.Add_Option_Argument("-3d", "force 3d mode");
    parse_args.Add_Option_Argument("-s", "write out slices as rgb files");
    parse_args.Add_Option_Argument("-b", "if -s selected, only write out binary slices");
    parse_args.Add_Option_Argument("-x", "along x direction");
    parse_args.Add_Option_Argument("-y", "along y direction");
    parse_args.Add_Option_Argument("-z", "along z direction");
    parse_args.Add_Option_Argument("-v", "display verbose information");
    parse_args.Set_Extra_Arguments(1, "<filename>");

    int extraarg = parse_args.Parse(argc, argv);

    if (extraarg < argc)
        strcpy(filename, argv[extraarg]);
    else
        return -1;

    // By default detection dimension from file extension
    if (Is_Phi_File(filename)) dimension = 3;
    else if (Is_Phi2D_File(filename)) dimension = 2;

    int slice_direction = 1;
    if (parse_args.Get_Option_Value("-x")) slice_direction = 1;
    if (parse_args.Get_Option_Value("-y")) slice_direction = 2;
    if (parse_args.Get_Option_Value("-z")) slice_direction = 3;

    if (parse_args.Get_Option_Value("-double")) type_double = true;
    if (parse_args.Get_Option_Value("-float")) type_double = false;
    if (parse_args.Get_Option_Value("-3d")) dimension = 3;
    if (parse_args.Get_Option_Value("-2d")) dimension = 2;

    if (dimension == -1) {
        std::cerr << "Could not determine dimension for file '" << filename << "'" << std::endl;
        return 1;
    }

    cout << "Filename: " << filename << " [" << ((type_double)?"double":"float") << ", " << dimension << "D]" << endl;

    if (dimension == 3)
    {
        if(!type_double) Print_Info_3D<float>(filename,
                                 parse_args.Get_Option_Value("-s"),
                                 parse_args.Get_Option_Value("-b"),
                                 slice_direction, parse_args.Get_Option_Value("-v"));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        else Print_Info_3D<double>(filename,
                                  parse_args.Get_Option_Value("-s"),
                                  parse_args.Get_Option_Value("-b"),
                                  slice_direction, parse_args.Get_Option_Value("-v"));
#else
        else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
    }
    else if (dimension == 2)
    {
        if(!type_double) Print_Info_2D<float>(filename, parse_args.Get_Option_Value("-v"));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        else Print_Info_2D<double>(filename, parse_args.Get_Option_Value("-v"));
#else
        else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
    }
}
