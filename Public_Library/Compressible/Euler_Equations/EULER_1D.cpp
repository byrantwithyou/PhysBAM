//#####################################################################
// Copyright 2002, 2003, Doug Enright, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Compressible/Euler_Equations/EULER.h>
#include <Compressible/Euler_Equations/EULER_1D.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void EULER_1D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.m;
    
    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_ghost(-2,m+3);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time);
    
    if(cut_out_grid) conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem_F);   
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,1> > psi(0,m);psi.Fill(true);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem_F);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T EULER_1D<T>::
CFL()
{
    int m=grid.m;T dx=grid.dx;
    
    ARRAY<T,VECTOR<int,1> > u_minus_c(0,m),u_plus_c(0,m);
    for(int i=0;i<m;i++){
        if(!cut_out_grid || (cut_out_grid && (*psi_pointer)(i))){
            T u=U(i)(1)/U(i)(0);
            T sound_speed=eos->c(U(i)(0),e(U(i)(0),U(i)(1),U(i)(2)));
            u_minus_c(i)=u-sound_speed;u_plus_c(i)=u+sound_speed;}}
    T dt_convect=max(u_minus_c.Max_Abs(),u_plus_c.Max_Abs())/dx;
    dt_convect=max(dt_convect,1/max_time_step);
    return 1/dt_convect;
}
//#####################################################################
}
