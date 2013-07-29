//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BURGERS_1D  
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Compressible/Burgers_Equation/BURGERS_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void BURGERS_1D<T>::
Euler_Step(const T dt,const T time)
{
    int m=grid.counts.x;
    int ghost_cells=3;

    ARRAY<TV,VECTOR<int,1> > U_ghost(-ghost_cells,m+ghost_cells);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);

    // doesn't  support cut out grids
    ARRAY<bool,VECTOR<int,1> > psi(0,m);psi.Fill(true);
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,TV>*,1> eigensystem(&eigensystem_F);
    conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T BURGERS_1D<T>::
CFL()
{
    int m=grid.counts.x;T dx=grid.dX.x;

    ARRAY<T,VECTOR<int,1> > u(0,m);
    for(int i=0;i<m;i++) u(i)=U(i)(0);
    T dt_convect=u.Max_Abs()/dx;

    return 1/max(dt_convect,1/max_time_step);
}
//#####################################################################
#if 0 // broken
template class BURGERS_1D<float>;
template class BURGERS_1D<double>;
#endif
