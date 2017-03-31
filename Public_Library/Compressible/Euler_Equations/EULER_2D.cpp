//#####################################################################
// Copyright 2002, 2003, Doug Enright, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Compressible/Euler_Equations/EULER.h>
#include <Compressible/Euler_Equations/EULER_2D.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void EULER_2D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.m,n=grid.n;
    
    ARRAY<TV_DIMENSION,VECTOR<int,2> > U_ghost(-2,m+3,-2,n+3);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time);
    
    if(cut_out_grid) conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem_F,eigensystem_G);  
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,2> > psi(0,m,0,n);psi.Fill(true);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem_F,eigensystem_G);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T EULER_2D<T>::
CFL()
{
    int m=grid.m,n=grid.n;T dx=grid.dx,dy=grid.dy;
    
    ARRAY<T,VECTOR<int,2> > u_minus_c(0,m,0,n),u_plus_c(0,m,0,n),v_minus_c(0,m,0,n),v_plus_c(0,m,0,n);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        if(!cut_out_grid || (cut_out_grid && (*psi_pointer)(i,j))){
            T u=U(1,i,j)/U(0,i,j),v=U(2,i,j)/U(0,i,j);
            T sound_speed=eos->c(U(0,i,j),e(U(0,i,j),U(1,i,j),U(2,i,j),U(3,i,j)));
            u_minus_c(i,j)=u-sound_speed;u_plus_c(i,j)=u+sound_speed;
            v_minus_c(i,j)=v-sound_speed;v_plus_c(i,j)=v+sound_speed;}}
    T dt_convect=max(u_minus_c.Max_Abs(),u_plus_c.Max_Abs())/dx+max(v_minus_c.Max_Abs(),v_plus_c.Max_Abs())/dy;
    dt_convect=max(dt_convect,1/max_time_step);
    return 1/dt_convect;
}           
}
