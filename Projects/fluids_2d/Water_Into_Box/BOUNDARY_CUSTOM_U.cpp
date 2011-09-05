//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Fedkiw - August 8, 2001
//#####################################################################
#include "BOUNDARY_CUSTOM_U.h"
using namespace PhysBAM;

//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
void BOUNDARY_CUSTOM_U::
Fill_Ghost_Cells(GRID_2D& grid,ARRAY_2D& u,ARRAY_2D& u_ghost,const double dt,const double time,const int number_of_ghost_cells)
{
    if(u.length != 1){Default();return;}

    int m=grid.m,n=grid.n;
    int i,j;

    put(u,u_ghost); // interior

    for(j=1;j<=n;j++){
        u_ghost(0,j)=-u_ghost(2,j);  // left 
        u_ghost(-1,j)=-u_ghost(3,j); // left 
        u_ghost(-2,j)=-u_ghost(4,j); // left 
        u_ghost(m+1,j)=-u_ghost(m-1,j);  // right 
        u_ghost(m+2,j)=-u_ghost(m-2,j);  // right 
        u_ghost(m+3,j)=-u_ghost(m-3,j);} // right 
    for(i=-2;i<=m+3;i++){
        u_ghost(i,0)=u_ghost(i,1);  // bottom 
        u_ghost(i,-1)=u_ghost(i,2); // bottom 
        u_ghost(i,-2)=u_ghost(i,3);} // bottom 
    // new - constant extrapolation on top
    for(i=-2;i<=m+3;i++){
        u_ghost(i,n+3)=u_ghost(i,n+2)=u_ghost(i,n+1)=u_ghost(i,n);} // top 
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
void BOUNDARY_CUSTOM_U::
Apply_Boundary_Condition(GRID_2D& grid,ARRAY_2D& u,const double time) 
{
    if(u.length != 1){Default();return;}
    
    int m=grid.m,n=grid.n;
    int j;

    for(j=1;j<=n;j++) u(1,j)=u(m,j)=0;  
} 
//#####################################################################

