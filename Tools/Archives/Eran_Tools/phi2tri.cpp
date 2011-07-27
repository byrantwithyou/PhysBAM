#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include "../../Public_Library/Geometry/LEVELSET_IMPLICIT_SURFACE.h"
#include <stdio.h>
#include <string.h>

using namespace std;
using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class T> void Phi_To_Tri(const char *input_filename,
                                  const char *output_filename,
                                  T resample_factor=1)
{
    LEVELSET_IMPLICIT_SURFACE<T> *implicit_surface = 0;
    FILE_UTILITIES::Create_From_File<T>(input_filename,implicit_surface);

    TRIANGULATED_SURFACE<T> *triangulated_surface = TRIANGULATED_SURFACE<T>::Create();

    if (resample_factor != 1)
    {
        GRID_3D<T> grid = implicit_surface->levelset.grid;
        T MAC_offset = grid.MAC_offset;
        grid.Initialize((int)(resample_factor*grid.m),(int)(resample_factor*grid.n),(int)(resample_factor*grid.mn),grid.xmin,grid.xmax,grid.ymin,grid.ymax,grid.zmin,grid.zmax);
        if (MAC_offset == 0.5) grid.Set_MAC_Grid();
        implicit_surface->levelset.Calculate_Triangulated_Surface_From_Marching_Tetrahedra(grid, *triangulated_surface);
    }
    else
    {
        implicit_surface->levelset.Calculate_Triangulated_Surface_From_Marching_Tetrahedra(*triangulated_surface);
    }

    ofstream out(output_filename, ios::binary);
    triangulated_surface->template Write<T>(out);

    delete implicit_surface;
    delete triangulated_surface;
}

int main(int argc, char *argv[])
{
    bool type_double = false;
    char filename[256];
    char output_filename[256];

    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Double_Argument("-resample", 1, "resample grid");
    parse_args.Add_String_Argument("-o", "", "output", "output file name");
    parse_args.Set_Extra_Arguments(1, "<filename>");

    int extraarg = parse_args.Parse(argc, argv);

    if (extraarg < argc)
        strcpy(filename, argv[extraarg]);
    else
        return -1;

    if (parse_args.Get_Option_Value("-float"))
        type_double = false;
    else if (parse_args.Get_Option_Value("-double"))
        type_double = true;

    if (!Is_Phi_File(filename))
    {
        cerr << "Not a phi file: " << filename << endl;
        return -1;
    }

    if (!parse_args.Is_Value_Set("-o"))
    {
        std::string basename=FILE_UTILITIES::Get_Basename(filename);
        sprintf(output_filename, "%s.tri", basename.c_str());
    }
    else
        strcpy(output_filename, parse_args.Get_String_Value("-o").c_str());

    cout << "Input filename: " << filename << endl;
    cout << "Output filename: " << output_filename << endl;

    if(!type_double) Phi_To_Tri<float>(filename, output_filename, parse_args.Get_Double_Value("-resample"));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else Phi_To_Tri<double>(filename, output_filename, parse_args.Get_Double_Value("-resample"));
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
}
