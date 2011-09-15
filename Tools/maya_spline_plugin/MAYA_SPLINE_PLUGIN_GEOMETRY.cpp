//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <iostream>
#include "MAYA_SPLINE_PLUGIN_GEOMETRY.h"

using namespace PhysBAM;

//#####################################################################
// Function MAYA_SPLINE_PLUGIN_GEOMETRY
//#####################################################################
MAYA_SPLINE_PLUGIN_GEOMETRY::
MAYA_SPLINE_PLUGIN_GEOMETRY()
{
    Initialize(GRID_3D<float>(2,2,2,0,1,0,1,0,1));
}
//#####################################################################
// Function ~MAYA_SPLINE_PLUGIN_GEOMETRY
//#####################################################################
MAYA_SPLINE_PLUGIN_GEOMETRY::
~MAYA_SPLINE_PLUGIN_GEOMETRY() 
{}
//#####################################################################
// Function Initialize
//#####################################################################
void MAYA_SPLINE_PLUGIN_GEOMETRY::
Initialize(const GRID_3D<float>& grid_input,const bool zeroed)
{ 
    grid=grid_input;
    controls.Resize(grid,1,false,false);
    if(zeroed) for(int i=0;i<=grid.m+1;i++) for(int j=0;j<=grid.n+1;j++) for(int ij=0;ij<=grid.mn+1;ij++) controls(i,j,ij)=VECTOR_3D<float>();
    else for(int i=0;i<=grid.m+1;i++) for(int j=0;j<=grid.n+1;j++) for(int ij=0;ij<=grid.mn+1;ij++){
        VECTOR_3D<float> value=grid.X(i,j,ij);controls(i,j,ij)=value;}
}
//#####################################################################
// Function Get_Number_Vertices
//#####################################################################
int MAYA_SPLINE_PLUGIN_GEOMETRY::
Get_Number_Vertices() 
{
    return (grid.m+2)*(grid.n+2)*(grid.mn+2);
}
//#####################################################################
