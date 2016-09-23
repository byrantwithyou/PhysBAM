//#####################################################################
// Copyright 2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CUSTOM
//#####################################################################
#ifndef __BOUNDARY_CUSTOM__
#define __BOUNDARY_CUSTOM__

#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class T,class TV>
class BOUNDARY_CUSTOM:public BOUNDARY<T,TV>
{
public:
    BOUNDARY_CUSTOM() 
    {}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<TV,VECTOR<int,1> >& u,ARRAY<TV,VECTOR<int,1> >& u_ghost,const T dt=0,const T time=0,const int number_of_ghost_cells=3) const {};
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<TV,VECTOR<int,2> >& u,ARRAY<TV,VECTOR<int,2> >& u_ghost,const T dt=0,const T time=0,const int number_of_ghost_cells=3) const ;
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<TV,VECTOR<int,3> >& u,ARRAY<TV,VECTOR<int,3> >& u_ghost,const T dt=0,const T time=0,const int number_of_ghost_cells=3) const {};
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<TV,VECTOR<int,1> >& u,const T time=0) const {};
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<TV,VECTOR<int,2> >& u,const T time=0) const ;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<TV,VECTOR<int,3> >& u,const T time=0) const {}; 
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T,class TV> void BOUNDARY_CUSTOM<T,TV>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<TV,VECTOR<int,2> >& u,ARRAY<TV,VECTOR<int,2> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    if(u.length != 1){assert(false);return;}

    int m=grid.m,n=grid.n;
    int i,j;

    ARRAY<TV,VECTOR<int,2> >::put(u,u_ghost); // interior

    for(j=0;j<n;j++){
        u_ghost(0,j)=u_ghost(m-1,j);  // left 
        u_ghost(-1,j)=u_ghost(m-2,j); // left 
        u_ghost(-2,j)=u_ghost(m-3,j); // left 
        u_ghost(m+1,j)=u_ghost(2,j);  // right 
        u_ghost(m+2,j)=u_ghost(3,j);  // right 
        u_ghost(m+3,j)=u_ghost(4,j);} // right 
    for(i=-2;i<=m+3;i++){
        u_ghost(i,0)=u_ghost(i,1);  // bottom 
        u_ghost(i,-1)=u_ghost(i,1); // bottom 
        u_ghost(i,-2)=u_ghost(i,1); // bottom 
        u_ghost(i,n+1)=u_ghost(i,n);  // top 
        u_ghost(i,n+2)=u_ghost(i,n);  // top 
        u_ghost(i,n+3)=u_ghost(i,n);} // top 
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T,class TV> void BOUNDARY_CUSTOM<T,TV>::
Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<TV,VECTOR<int,2> >& u,const T time) const
{
    if(u.length != 1){assert(false);return;}
    
    int m=grid.m,n=grid.n;
    int j;

    for(j=0;j<n;j++) u(1,j)=u(m,j);  // left 
    //for(i=0;i<m;i++) u(i,1)=0;  // top 
}
//#####################################################################
}
#endif
