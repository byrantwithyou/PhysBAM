//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/max.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Compressible/Reactive_Euler_Equations/REACTIVE_EULER_2D.h>
namespace PhysBAM{
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T> void REACTIVE_EULER_2D<T>::
Euler_Step(const T dt,const T time)
{  
    int ghost_cells=3;
    
    ARRAY<TV_DIMENSION,VECTOR<int,2> > U_ghost(grid.Domain_Indices(ghost_cells));
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);

    ARRAY<bool,FACE_INDEX<TV::m> > psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,5>*,2> eigensystem(&eigensystem_F,&eigensystem_G);
    if(cut_out_grid) 
        conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem,eigensystem,psi_N,face_velocities);
    else{ // not a cut out grid
        ARRAY<bool,VECTOR<int,2> > psi(grid.counts);
        psi.Fill(1);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T REACTIVE_EULER_2D<T>::
CFL()
{
    int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    
    ARRAY<T,VECTOR<int,2> > u_minus_c(grid.counts),u_plus_c(grid.counts),v_minus_c(grid.counts),v_plus_c(grid.counts);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        if(!cut_out_grid || (cut_out_grid && (*psi_pointer)(i,j)==1)){
            T u=U(i,j)(1)/U(i,j)(0),v=U(i,j)(2)/U(i,j)(0);
            T Y=U(i,j)(4)/U(i,j)(0);
            T sound_speed=eos.c(U(i,j)(0),e(U(i,j)(0),U(i,j)(1),U(i,j)(2),U(i,j)(3)),Y);
            u_minus_c(i,j)=u-sound_speed;u_plus_c(i,j)=u+sound_speed;
            v_minus_c(i,j)=v-sound_speed;v_plus_c(i,j)=v+sound_speed;}}
    T dt_convect=max(u_minus_c.Max_Abs(),u_plus_c.Max_Abs())/dx+max(v_minus_c.Max_Abs(),v_plus_c.Max_Abs())/dy;
    return 1/dt_convect;
}
//#####################################################################
template class REACTIVE_EULER_2D<float>;
template class REACTIVE_EULER_2D<double>;
}
