//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Spherical_Euler_Equations/SPHERICAL.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void SPHERICAL<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x;
    int i,k;
    int ghost_cells=3;

    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_ghost(-ghost_cells,m+ghost_cells);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);
    
    // evaluate source terms
    ARRAY<TV_DIMENSION,VECTOR<int,1> > S(0,m);
    for(i=0;i<m;i++){
        T rho=U(i)(0);
        T u=U(i)(1)/U(i)(0);
        T e=U(i)(2)/U(i)(0)-sqr(u)/2;
        T rho_u=U(i)(1);
        T coefficient=-2/grid.X(VECTOR<int,1>(i)).x;
        S(i)(0)=coefficient*rho_u;
        S(i)(1)=coefficient*rho_u*u;
        S(i)(2)=coefficient*(U(i)(2)+eos->p(rho,e))*u;}
    
    T_FACE_ARRAYS_BOOL psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    T_FACE_ARRAYS_SCALAR face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,3> >*,1> eigensystem(&eigensystem_F);
    if(cut_out_grid) conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem,eigensystem,psi_N,face_velocities);
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,1> > psi(0,m);psi.Fill(1);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);}

    // add source terms
    for(i=0;i<m;i++) for(k=0;k<3;k++) U(i)(k)+=dt*S(i)(k);
    
    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
template class SPHERICAL<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SPHERICAL<double>;
#endif
