//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
using namespace PhysBAM;

typedef float T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> TV_INT;
typedef VECTOR<T,2> V2;

int main(int argc, char* argv[])
{
    std::string file="surface.tri",input="<none>";
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-o",&file,"file","output filename");
    parse_args.Add("-i",&input,"file","input level set");
    parse_args.Parse();

    TRIANGULATED_SURFACE<T>& ts=*TRIANGULATED_SURFACE<T>::Create();
    if(input=="<none>"){
        RANDOM_NUMBERS<T> random;
        GRID<TV> grid(TV_INT()+50,RANGE<TV>::Centered_Box());
        ARRAY<T,TV_INT> phi(grid.Node_Indices());
        phi.Fill(1);
        for(NODE_ITERATOR<TV> it(grid,-1);it.Valid();it.Next())
            phi(it.index)=random.Get_Uniform_Number(-1,1);
        MARCHING_CUBES<TV>::Create_Surface(ts,grid,phi);}
    else{
        LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* io;
        FILE_UTILITIES::Create_From_File<T>(input,io);
        ARRAY<T,TV_INT>& phi=io->levelset.phi;
        MARCHING_CUBES<TV>::Create_Surface(ts,io->levelset.grid,phi);
        delete io;}

    ts.mesh.Initialize_Boundary_Mesh();
    LOG::cout<<ts.mesh.boundary_mesh->elements.m<<std::endl;
    FILE_UTILITIES::Write_To_File<T>(file,ts);

    delete &ts;
    return 0;
}
