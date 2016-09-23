#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <cstring>
#include <fstream>
using namespace PhysBAM;
using namespace std;
//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    typedef float T;
    typedef VECTOR<T,3> TV;

    typedef float RW;
    STREAM_TYPE stream_type((RW()));

    if(argc<3){
        LOG::cout<<"Usage ./phi2tri <phi input> <tri output>"<<std::endl;
        return 1;}
    
    GRID<VECTOR<T,3> > grid;
    ARRAY<T,VECTOR<int,3> > phi;
    LEVELSET_IMPLICIT_OBJECT<TV> implicit_surface(grid,phi);
    FILE_UTILITIES::Read_From_File(stream_type,argv[1],implicit_surface);
    std::cout<<"Grid : "<<grid<<std::endl;

    TRIANGULATED_SURFACE<T>& triangulated_surface=*DUALCONTOUR_3D<T>::Create_Triangulated_Surface_From_Levelset(implicit_surface.levelset);

    std::cout<<"Surface triangles : "<<triangulated_surface.mesh.elements.m<<std::endl;
    std::cout<<"Surface particles : "<<triangulated_surface.particles.Size()<<std::endl;

    FILE_UTILITIES::Write_To_File(stream_type,argv[2],triangulated_surface);

    return 0;
}
//#################################################################
