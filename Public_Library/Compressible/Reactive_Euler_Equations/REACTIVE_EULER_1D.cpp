//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/maxabs.h>
#include <Core/Math_Tools/sqr.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Compressible/Reactive_Euler_Equations/REACTIVE_EULER_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void REACTIVE_EULER_1D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x;
    int ghost_cells=3;
    
    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_ghost(-ghost_cells,m+ghost_cells);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);
    
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,1> eigensystem(&eigensystem_F);
    if(cut_out_grid) conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem,eigensystem,psi_N,face_velocities);
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,1> > psi(0,m);psi.Fill(1);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T REACTIVE_EULER_1D<T>::
CFL()
{
    int m=grid.counts.x;T dx=grid.dX.x;
    
    ARRAY<T,VECTOR<int,1> > u_minus_c(0,m),u_plus_c(0,m);
    for(int i=0;i<m;i++){
        if(!cut_out_grid || (cut_out_grid && (*psi_pointer)(i)==1)){
            T u=U(i)(1)/U(i)(0);
            T Y=U(i)(3)/U(i)(0);
            T sound_speed=eos.c(U(i)(0),e(U(i)(0),U(i)(1),U(i)(2)),Y);
            u_minus_c(i)=u-sound_speed;u_plus_c(i)=u+sound_speed;}}
    T dt_convect=max(u_minus_c.Max_Abs(),u_plus_c.Max_Abs())/dx;
    return 1/dt_convect;
}
//#####################################################################
#if 0 // broken
template class REACTIVE_EULER_1D<float>;
template class REACTIVE_EULER_1D<double>;
#endif
