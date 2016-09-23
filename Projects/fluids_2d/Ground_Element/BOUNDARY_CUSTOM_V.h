//#####################################################################
// Copyright 2002, Ronald Fedkiw, Duc Nguyen
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CUSTOM_V
//#####################################################################
//
//#####################################################################
// Fedkiw - August 8, 2001
// Nguyen - August 30, 2002
//#####################################################################
#ifndef __BOUNDARY_CUSTOM_V__
#define __BOUNDARY_CUSTOM_V__

#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class T,class T2>
class BOUNDARY_CUSTOM_V:public BOUNDARY<T,T2>
{
public:
    BOUNDARY_CUSTOM_V() 
                                         :BOUNDARY<T,T2>()
    {}

//#####################################################################
    void Fill_Ghost_Cells(GRID_2D<T>& grid,ARRAY<T2,VECTOR<int,2> >& u,ARRAY<T2,VECTOR<int,2> >& u_ghost,const T dt=0,const T time=0,const int number_of_ghost_cells=3) const;
    void Apply_Boundary_Condition(GRID_2D<T>& grid,ARRAY<T2,VECTOR<int,2> >& u,const T time=0) const;
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T,class T2> void BOUNDARY_CUSTOM_V<T,T2>::
Fill_Ghost_Cells(GRID_2D<T>& grid,ARRAY<T2,VECTOR<int,2> >& u,ARRAY<T2,VECTOR<int,2> >& u_ghost,const T time) const
{
    if(u.length != 1){Default();return;}

    int m=grid.m,n=grid.n;
    int i,j;

    ARRAY<T2,VECTOR<int,2> >::put(u,u_ghost); // interior

    // constant boundary condition for left and right
    for(j=0;j<n;j++){
        u_ghost(-2,j)=u_ghost(-1,j)=u_ghost(0,j)=u_ghost(1,j);  // left 
        u_ghost(m+3,j)=u_ghost(m+2,j)=u_ghost(m+1,j)=u_ghost(m,j);} // right 
    // constant extrapolation on top
    for(i=-2;i<=m+3;i++){
        u_ghost(i,n+3)=u_ghost(i,n+2)=u_ghost(i,n+1)=u_ghost(i,n);} // top 
    for(i=-2;i<=m+3;i++){
        u_ghost(i,0)=-u_ghost(i,2);  // bottom 
        u_ghost(i,-1)=-u_ghost(i,3); // bottom 
        u_ghost(i,-2)=-u_ghost(i,4);} // bottom 
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T,class T2> void BOUNDARY_CUSTOM_V<T,T2>::
Apply_Boundary_Condition(GRID_2D<T>& grid,ARRAY<T2,VECTOR<int,2> >& u,const T time) const
{
    if(u.length != 1){Default();return;}
    
    int m=grid.m,n=grid.n;
    int i;

    for(i=0;i<m;i++) u(i,1)=(T2)0; // new - do not set top to zero u(i,n,ij)=0; 
}  
//#####################################################################
}
#endif
