#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/UNION_FIND.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Images/EPS_FILE.h>
#include <Geometry/Images/TEX_FILE.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include "MARCHING_TETRAHEDRA_CUTTING.h"
using namespace PhysBAM;

#define P(x) LOG::cout<<#x<<" : "<<(x)<<std::endl;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef float RW;
    typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<int,4> E;
    typedef VECTOR<T,3> TV;
    typedef MARCHING_TETRAHEDRA_CUTTING<TV>::DATA DATA;

    T w=(T).2,d=(T).4;
    TV_INT size(TV_INT()+4);
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-n",&size,"value","Domain size");
    parse_args.Add("-w",&w,"value","Slit width");
    parse_args.Add("-d",&d,"value","Slit depth");
    parse_args.Parse();

    GRID<TV> grid(size+1,RANGE<TV>(-TV(size)/2,TV(size)/2));
    TETRAHEDRALIZED_VOLUME<T> tv;
    tv.Initialize_Cube_Mesh_And_Particles(grid);

    ARRAY<T> phi0(tv.particles.X.m),phi1(tv.particles.X.m);
    phi0.Fill(-1);
    for(int i=0;i<tv.particles.X.m;i++){
        TV Z=tv.particles.X(i)-TV(0,0,size(2)/2);
        phi1(i)=sqr(Z.x/(w*size.x))+sqr(Z.z/(d*size.z))-1;}

    ARRAY<E> out_mesh[2];
    ARRAY<DATA> data[2];
    ARRAY<TV_INT> surface[2];
    ARRAY<int> node_map;
    ARRAY<TV> X(tv.particles.X);
    MARCHING_TETRAHEDRA_CUTTING<TV>::Fracture_Cutting(tv.mesh.elements,X,phi0,phi1,out_mesh,data,surface,node_map);
    const char* files[2]={"in.tri.gz","out.tri.gz"};
    for(int s=0;s<2;s++){
        TRIANGULATED_SURFACE<T> ts;
        ts.particles.Add_Elements(X.m);
        ts.particles.X=X;
        ts.mesh.elements=surface[s];
        ts.Update_Number_Nodes();
        FILE_UTILITIES::Write_To_File<RW>(files[s],ts);}
    FILE_UTILITIES::Write_To_File<RW>("notch.tet.gz",tv);

    return 0;
}
