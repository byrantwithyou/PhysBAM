//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL
//#####################################################################
//
// Assumes a MAC grid straddling the centerline on the left, applying reflection.
// Extrapolation is applied to the right, top and bottom.
//
//#####################################################################
#ifndef __BOUNDARY_EULER_EQUATIONS_CYLINDRICAL__
#define __BOUNDARY_EULER_EQUATIONS_CYLINDRICAL__

#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class T_input>
class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL:public BOUNDARY<VECTOR<T_input,2>,VECTOR<T_input,4> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<T,4> TV_DIMENSION;typedef VECTOR<int,2> TV_INT;
    enum {d=4};
public:
    BOUNDARY_EULER_EQUATIONS_CYLINDRICAL()
    {}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,2> >& u,ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,2> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,2> >& u,const T time) const override {} // do nothing
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T> void BOUNDARY_EULER_EQUATIONS_CYLINDRICAL<T>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,2> >& u,ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,2> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    int m=grid.counts.x,n=grid.counts.y;
    int i,j;

    ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,2> >::Put(u,u_ghost); // interior

    // left
    for(i=-3;i<0;i++) for(j=0;j<n;j++){
        T rho=u_ghost(-i,j)(0);
        T u_velocity=-u_ghost(-i,j)(1)/u_ghost(-i,j)(0);
        T v_velocity=u_ghost(-i,j)(2)/u_ghost(-i,j)(0);
        T e=u_ghost(-i,j)(3)/u_ghost(-i,j)(0)-(sqr(u_velocity)+sqr(v_velocity))/2;
        u_ghost(i,j)(0)=rho;
        u_ghost(i,j)(1)=rho*u_velocity;
        u_ghost(i,j)(2)=rho*v_velocity;
        u_ghost(i,j)(3)=rho*(e+(sqr(u_velocity)+sqr(v_velocity))/2);}

    // constant extrapolation
    for(j=0;j<n;j++) u_ghost(m+3,j)=u_ghost(m+2,j)=u_ghost(m+1,j)=u_ghost(m,j); // right
    for(i=-3;i<m+3;i++){
        u_ghost(i,-2)=u_ghost(i,-1)=u_ghost(i,0)=u_ghost(i,1);           // bottom
        u_ghost(i,n+3)=u_ghost(i,n+2)=u_ghost(i,n+1)=u_ghost(i,n);} // top
}
//#####################################################################
}
#endif
