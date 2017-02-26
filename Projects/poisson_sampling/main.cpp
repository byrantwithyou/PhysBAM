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
    vd.Init(RANGE<TV>::Unit_Box(),(T).001);
    vd.Visualize_State("Initial");

    T sample_density=10000*0;
    int n=ceil(sample_density*((vd.pieces.m?vd.pieces(0).subtree_area:0)+(vd.clipped_pieces.m?vd.clipped_pieces(0).subtree_area:0)));
    for(int i=0;i<n;i++){
        int p=vd.Choose_Piece();
        if(p<0) break;
        TV X=vd.Choose_Feasible_Point(p);
        Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));}
    vd.Visualize_State("Initial Sample");
    
    for(int i=0;;i++){
        T area=(vd.pieces.m?vd.pieces(0).subtree_area:0)+(vd.clipped_pieces.m?vd.clipped_pieces(0).subtree_area:0);
        printf("%i %g\n",i,area);
        int p=vd.Choose_Piece();
        if(p<0) break;
        TV X=vd.Choose_Feasible_Point(p);
        vd.Insert_Point(p,X);

        int n=ceil(sample_density*((vd.pieces.m?vd.pieces(0).subtree_area:0)+(vd.clipped_pieces.m?vd.clipped_pieces(0).subtree_area:0)));
        for(int i=0;i<n;i++){
            int p=vd.Choose_Piece();
            if(p<0) break;
            TV X=vd.Choose_Feasible_Point(p);
            Add_Debug_Particle(X,VECTOR<T,3>(0,1,0));}}
    vd.Visualize_State("End");
    vd.Sanity_Checks();
    
    Flush_Frame<TV>("end");

    return 0;
}

