//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Math_Tools/sqr.h>
#include <Compressible/Euler_Equations/CYLINDRICAL.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void CYLINDRICAL<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x,n=grid.counts.y;
    int ghost_cells=3;

    ARRAY<TV_DIMENSION,VECTOR<int,2> > U_ghost(grid.Domain_Indices(ghost_cells));
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);

    // evaluate source terms
    ARRAY<TV_DIMENSION,VECTOR<int,2> > S(grid.counts);
    for(RANGE_ITERATOR<2> it(S.domain);it.Valid();it.Next()){
        T rho=U(it.index)(0);
        T u=U(it.index)(0)/U(it.index)(0);
        T v=U(it.index)(1)/U(it.index)(0);
        T e=U(it.index)(2)/U(it.index)(0)-(sqr(u)+sqr(v))/2;
        T rho_u=U(it.index)(1);
        T coefficient=-1/grid.Node(it.index).x;
        S(it.index)(0)=coefficient*rho_u;
        S(it.index)(1)=coefficient*rho_u*u;
        S(it.index)(2)=coefficient*rho_u*v;
        S(it.index)(3)=coefficient*(U(it.index)(3)+eos->p(rho,e))*u;}
    
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,2> eigensystem(&eigensystem_F,&eigensystem_G);
    if(cut_out_grid) conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem,eigensystem,psi_N,face_velocities);
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,2> > psi(0,m,0,n);psi.Fill(1);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);}

    // add source terms
    for(RANGE_ITERATOR<2> it(S.domain);it.Valid();it.Next()) U(it.index)+=dt*S(it.index);
    
    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
namespace PhysBAM{
template class CYLINDRICAL<float>;
template class CYLINDRICAL<double>;
}
