//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
using namespace PhysBAM;

typedef float T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> TV_INT;
typedef VECTOR<T,2> V2;

int main(int argc, char* argv[])
{
    Initialize_Particles();
    Initialize_Read_Write_General_Structures();
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-o","surface.tri","output filename");
    parse_args.Add_String_Argument("-i","<none>","input level set");
    parse_args.Parse(argc,argv);
    std::string file=parse_args.Get_String_Value("-o");
    std::string input=parse_args.Get_String_Value("-i");

    TRIANGULATED_SURFACE<T>& ts=*TRIANGULATED_SURFACE<T>::Create();
    if(input=="<none>"){
        RANDOM_NUMBERS<T> random;
        GRID<TV> grid(TV_INT()+50,RANGE<TV>::Centered_Box());
        ARRAY<T,TV_INT> phi(grid.Node_Indices());
        phi.Fill(1);
        for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,-1);it.Valid();it.Next())
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
