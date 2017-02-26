//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include "VORONOI_DIAGRAM.h"

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,TV::m> TV_INT;

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    T radius=1./200;
    GRID<TV> grid(TV_INT()+1,RANGE<TV>::Unit_Box(),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),grid,"output");

    ARRAY<TV> X,Y;

    auto c1=clock();
    VORONOI_DIAGRAM<T> vd;
    vd.Sample_Fully(RANGE<TV>::Unit_Box(),radius);
    vd.Get_Samples(X);
    auto c2=clock();
    POISSON_DISK<TV> poisson_disk(radius);
    RANDOM_NUMBERS<T> r;
    poisson_disk.Sample(r,RANGE<TV>::Unit_Box(),Y);
    auto c3=clock();
    LOG::printf("%P %P\n",(c2-c1)/(double)CLOCKS_PER_SEC,(c3-c2)/(double)CLOCKS_PER_SEC);
    LOG::printf("%P %P\n",X.m,Y.m);

    // for(int i=0;i<X.m;i++)
    //     Add_Debug_Particle(X(i),VECTOR<T,3>(1,0,0));
    // Flush_Frame<TV>("voronoi");

    // for(int i=0;i<Y.m;i++)
    //     Add_Debug_Particle(Y(i),VECTOR<T,3>(1,0,0));
    // Flush_Frame<TV>("poisson disc");
    
    // Flush_Frame<TV>("end");

    return 0;
}

