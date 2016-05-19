//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    GRID<TV> grid(TV_INT()+2,RANGE<TV>::Unit_Box(),false);
    ARRAY<T,TV_INT> phi(TV_INT()+2);

    for(int c=0;c<256;c++)
    {
        LOG::printf("CASE: %i\n",c);
        for(int i=0;i<8;i++)
        {
            TV_INT j(i&1,i/2&1,i/4&1);
            phi(j)=(c&(1<<i))?-1:1;
        }
        TETRAHEDRALIZED_VOLUME<T> tv;
        MARCHING_CUBES<TV>::Create_Interior(tv,grid,phi,true);
        for(int i=0;i<tv.mesh.elements.m;i++)
        {
            LOG::printf("%P\n",tv.particles.X.Subset(tv.mesh.elements(i)));
        }
    }

    return 0;
}

