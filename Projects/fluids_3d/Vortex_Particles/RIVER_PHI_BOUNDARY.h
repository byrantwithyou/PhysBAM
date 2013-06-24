//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIVER_PHI_BOUNDARY  
//#####################################################################
#ifndef __RIVER_PHI_BOUNDARY__
#define __RIVER_PHI_BOUNDARY__

#include <Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{
template<class T,class T2>
class RIVER_PHI_BOUNDARY:public BOUNDARY<T,T>
{
public:
    LEVELSET_IMPLICIT_SURFACE<T>* implicit;
    T inflow_height;
    
    RIVER_PHI_BOUNDARY()
        :implicit(0),inflow_height(0)
    {}
    
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<T2,VECTOR<int,3> >& u,ARRAY<T2,VECTOR<int,3> >& u_ghost,const T dt=0,const T time=0,const int number_of_ghost_cells=3) const
    {ARRAY<T2,VECTOR<int,3> >::put(u,u_ghost);Fill_Left_Ghost_Cells(grid,u_ghost,time);Fill_Right_Ghost_Cells(grid,u_ghost,time);Fill_Bottom_Ghost_Cells(grid,u_ghost,time);Fill_Top_Ghost_Cells(grid,u_ghost,time);
    Fill_Front_Ghost_Cells(grid,u_ghost,time);Fill_Back_Ghost_Cells(grid,u_ghost,time);
    const int n=grid.n,mn=grid.mn;int j,ij;
    // override at the inflow
    for(j=0;j<n;j++) for(ij=0;ij<mn;ij++){
        if((implicit && implicit->Lazy_Inside(grid.X(1,j,ij))) || grid.y(j) > inflow_height) u_ghost(-2,j,ij)=u_ghost(-1,j,ij)=u_ghost(0,j,ij)=grid.dx;
        else u_ghost(-2,j,ij)=u_ghost(-1,j,ij)=u_ghost(0,j,ij)=-grid.dx;}
    // don't suck water out of outflow
    for(j=0;j<n;j++) for(ij=0;ij<mn;ij++) u_ghost(grid.m+1,j,ij)=u_ghost(grid.m+2,j,ij)=u_ghost(grid.m+3,j,ij)=u_ghost(grid.m,j,ij);}
//#####################################################################
};
}
#endif
