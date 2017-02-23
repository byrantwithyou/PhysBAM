//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
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
/*
    vd.First_Three_Points(TV(0,0),TV(1,0),TV(.5,1));
    vd.Visualize_State("Initial");
    vd.Sanity_Checks();

    while(vd.pieces.m){
        int p=vd.Choose_Piece();
        TV X=vd.pieces(p).h.Choose_Feasible_Point(vd.random,vd.radius);
        vd.Insert_Point(vd.pieces(p).coedge,X);
        vd.Visualize_State("After insert");}
*/
    vd.Init(RANGE<TV>::Unit_Box(),(T).1);
    vd.Visualize_State("Initial");

    for(int i=0;i<1;i++){
        int p=vd.Choose_Piece();
        TV X=vd.Choose_Feasible_Point(p);
        vd.Insert_Point(p,X);
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));
        vd.Visualize_State("After insertion");}

    for(int i=0;i<10000;i++){
        int p=vd.Choose_Piece();
        TV X=vd.Choose_Feasible_Point(p);
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));}
    vd.Visualize_State("End");

    Flush_Frame<TV>("end");

    return 0;
}

