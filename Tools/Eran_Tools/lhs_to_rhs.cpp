#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include <string>
#include "../../Public_Library/Geometry/LEVELSET_IMPLICIT_SURFACE.h"

using namespace PhysBAM;

int main()
{
    LEVELSET_IMPLICIT_SURFACE<double> *canyon=0;
    FILE_UTILITIES::Create_From_File<double>("canyon.phi",canyon);

    GRID_3D<double> &grid = canyon->levelset.grid;
    grid.zmin *= -1; grid.zmax *= -1;
    double tmp = grid.zmin; grid.zmin = grid.zmax; grid.zmax = tmp;

    for (int i = 1; i <= grid.m; i++) for (int j = 1; j <= grid.n; j++)
    {
        for (int k = 1; k <= grid.mn/2; k++) 
        {
            tmp = canyon->levelset.phi(i,j,k);
            canyon->levelset.phi(i,j,k) = canyon->levelset.phi(i,j,grid.mn-k+1);
            canyon->levelset.phi(i,j,grid.mn-k+1) = tmp;
        }
    }

    std::ofstream output_file("canyon_flipped.phi", std::ios::binary);
    canyon->Write<double>(output_file);
    output_file.close();
}
