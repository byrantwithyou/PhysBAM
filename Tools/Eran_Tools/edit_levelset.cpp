#include <PhysBAM_Tools/Images/RGB_FILE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <fstream>
#include <iostream>
#include "../../Public_Library/Geometry/LEVELSET_IMPLICIT_SURFACE.h"
#include <stdio.h>
#include <string.h>

using namespace PhysBAM;
using namespace FILE_UTILITIES;

template<class T> void Do_It(const std::string& filename,PARSE_ARGS& parse_args)
{
    std::string reference_filename="levelset.244";
    LEVELSET_IMPLICIT_SURFACE<T> *reference_implicit_surface=0;
    FILE_UTILITIES::Create_From_File<T>(reference_filename,reference_implicit_surface);
    GRID_3D<T> &reference_grid = reference_implicit_surface->levelset.grid;
    ARRAYS<VECTOR<T,3> > &reference_phi = reference_implicit_surface->levelset.phi;

    T y_cutoff=8;

    LEVELSET_IMPLICIT_SURFACE<T> *implicit_surface=0;
    FILE_UTILITIES::Create_From_File<T>(filename,implicit_surface);
    GRID_3D<T> &grid = implicit_surface->levelset.grid;
    ARRAYS<VECTOR<T,3> > &phi = implicit_surface->levelset.phi;

#if 0
    T domain_scale=10;
    CYLINDER<T> cylinder(domain_scale*VECTOR<T,3>(.55,.9,.5),domain_scale*VECTOR<T,3>(.55,1,.5),domain_scale*.05);
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int k=0;k<grid.mn;k++)
        if(cylinder.Inside(grid.X(i,j,k),3*grid.min_dx_dy_dz)) phi(i,j,k)=min(phi(i,j,k),cylinder.Signed_Distance(grid.X(i,j,k)));
#endif
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int k=0;k<grid.mn;k++)
        if(grid.y(j)>8) phi(i,j,k)=reference_phi(i,j,k);

    FILE_UTILITIES::Write_To_File<T>(filename+".out",*implicit_surface);
}

int main(int argc, char *argv[])
{
    char filename[256];

    PARSE_ARGS parse_args;

    
    parse_args.Set_Extra_Arguments(1, "<filename>");
    int extraarg = parse_args.Parse(argc, argv);
    if (extraarg < argc)
        strcpy(filename, argv[extraarg]);
    else
        return -1;

    Do_It<float>(filename,parse_args);
}
