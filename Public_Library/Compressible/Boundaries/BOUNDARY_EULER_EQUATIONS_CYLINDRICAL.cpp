//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_CYLINDRICAL.h>
namespace PhysBAM{
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
template class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL<double>;
template class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL<float>;
}
