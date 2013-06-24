#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/UNION_FIND.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Images/EPS_FILE.h>
#include <Geometry/Images/TEX_FILE.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "MARCHING_TETRAHEDRA_CUTTING.h"
using namespace PhysBAM;

#define P(x) LOG::cout<<#x<<" : "<<(x)<<std::endl;

int main(int argc,char* argv[])
{
    typedef double T;
    typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<int,4> E;
    typedef VECTOR<T,3> TV;
    typedef MARCHING_TETRAHEDRA_CUTTING<TV>::DATA DATA;

    int seed=time(0),size=4;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-s",&seed,"value","Seed");
    parse_args.Add("-n",&size,"value","Domain size");
    parse_args.Parse();

    ARRAY<E> in_mesh;
    in_mesh.Append(E(1,2,5,0));
    in_mesh.Append(E(2,4,5,0));
    in_mesh.Append(E(4,3,5,0));
    in_mesh.Append(E(3,1,5,0));
    ARRAY<T> phi0;
    phi0.Append(1.5);
    phi0.Append(0.5);
    phi0.Append(0.5);
    phi0.Append(0.5);
    phi0.Append(0.5);
    phi0.Append(-0.5);
    ARRAY<T> phi1;
    phi1.Append(-1);
    phi1.Append(-1);
    phi1.Append(-1);
    phi1.Append(-1);
    phi1.Append(-1);
    phi1.Append(-1);
    ARRAY<TV> X;
    X.Append(TV(0,0,-1));
    X.Append(TV(-1,-1,0));
    X.Append(TV(1,-1,0));
    X.Append(TV(-1,1,0));
    X.Append(TV(1,1,0));
    X.Append(TV(0,0,1));

    ARRAY<E> out_mesh[2];
    ARRAY<DATA> data[2];
    ARRAY<TV_INT> surface[2];
    ARRAY<int> node_map;
    MARCHING_TETRAHEDRA_CUTTING<TV>::Fracture_Cutting(in_mesh,X,phi0,phi1,out_mesh,data,surface,node_map);
    P(in_mesh);
    P(X);
    P(phi0);
    P(phi1);
    P(out_mesh[0]);
    P(out_mesh[1]);
    P(surface[0]);
    P(surface[1]);
    P(node_map);
    for(int s=0;s<2;s++){
        P(s);
        for(int i=0;i<data[s].m;i++){
            P(i);
            P(data[s](i).parent);
            P(data[s](i).element);
            P(data[s](i).volume);}}

    return 0;
}
