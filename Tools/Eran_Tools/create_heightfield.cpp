#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace PhysBAM;

int main(int argc, char *argv[])
{
    double xmin, xmax, ymin, ymax;
    int m, n;
    int i, j;

    m = 64; n = 64;
    xmin = -859.209; xmax = 695.791; ymin = -187.55; ymax = 282.45;

    if (argc > 1) m = atoi(argv[1]);
    if (argc > 2) n = atoi(argv[2]);
    if (argc > 3) xmin = atof(argv[3]);
    if (argc > 4) xmax = atof(argv[4]);
    if (argc > 5) ymin = atof(argv[5]);
    if (argc > 6) ymax = atof(argv[6]);

    std::cout << "Doing " << m << " x " << n << ", " << xmin << " " << xmax << " " << ymin << " " << ymax << std::endl;

    GRID_2D<float> grid(m,n,xmin,xmax,ymin,ymax);
    ARRAYS<VECTOR<float,2> > ground(1,m,1,n);
    
    TRIANGULATED_SURFACE<double> *canyon=0;
    FILE_UTILITIES::Create_From_File<double>("canyon.tri",canyon);

    canyon->Update_Bounding_Box();
    canyon->Initialize_Triangle_Hierarchy();

    for (i=1;i<=m;i++) for (j=1;j<=n;j++)
    {
        std::cout << i << "," << j << " " << std::flush;
        RAY_3D<double> ray(VECTOR<double,3>(grid.x(i),1000,grid.y(j)), VECTOR<double,3>(0,-1,0));
        if (canyon->Intersection(ray)) {
            ground(i,j) = ray.Point(ray.t_max).y;
        }
        else
        {
            ground(i,j) = -155.936;
        }
    }

    std::ofstream output_file;
    output_file.open("canyon.heightfield", std::ios::binary); ground.Write<float>(output_file); output_file.close();
    output_file.open("canyon.grid", std::ios::binary); grid.Write<float>(output_file); output_file.close();
}
