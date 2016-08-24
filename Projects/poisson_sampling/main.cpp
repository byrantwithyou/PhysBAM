//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
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

    GRID<TV> grid(TV_INT()+1,RANGE<TV>::Unit_Box(),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),grid,"output");
    
    VORONOI_DIAGRAM<T> vd;

    vd.First_Three_Points(TV(0,0),TV(1,0),TV(.5,1));
    vd.Visualize_State("Initial");
    vd.Sanity_Checks();

    while(vd.pieces.m){
        int p=vd.Choose_Piece();
        TV X=vd.Choose_Feasible_Point(vd.pieces(p).coedge);
        vd.Insert_Point(vd.pieces(p).coedge,X);
        vd.Visualize_State("After insert");}

    Flush_Frame<TV>("end");

    return 0;
}

