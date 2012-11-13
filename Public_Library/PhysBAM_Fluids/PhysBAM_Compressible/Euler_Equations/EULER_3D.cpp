//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void EULER_3D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;
    int ghost_cells=3;
    
    ARRAY<TV_DIMENSION,TV_INT> U_ghost(-ghost_cells,m+ghost_cells,-ghost_cells,n+ghost_cells,1-ghost_cells,mn+ghost_cells);
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);
    
    T_FACE_ARRAYS_BOOL psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    T_FACE_ARRAYS_SCALAR face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,VECTOR<T,5> >*,3> eigensystem(&eigensystem_F,&eigensystem_G,&eigensystem_H);
    if(psi_pointer) 
        conservation->Update_Conservation_Law(grid,U,U_ghost,*psi_pointer,dt,eigensystem,eigensystem,psi_N,face_velocities);
    else{ // not a cut out grid
        ARRAY<bool,TV_INT> psi(grid.counts);
        psi.Fill(1);
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);}

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T EULER_3D<T>::
CFL()
{
    ARRAY<T,TV_INT> u_minus_c(grid.counts),u_plus_c(grid.counts),v_minus_c(grid.counts),
             v_plus_c(grid.counts),w_minus_c(grid.counts),w_plus_c(grid.counts);
    for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()){
        if(!psi_pointer || (*psi_pointer)(it.index)==1){
            T u=U(it.index)(1)/U(it.index)(0),v=U(it.index)(2)/U(it.index)(0),w=U(it.index)(3)/U(it.index)(0);
            T sound_speed=eos->c(U(it.index)(0),e(U(it.index)(0),U(it.index)(1),U(it.index)(2),U(it.index)(3),U(it.index)(4)));
            u_minus_c(it.index)=u-sound_speed;u_plus_c(it.index)=u+sound_speed;
            v_minus_c(it.index)=v-sound_speed;v_plus_c(it.index)=v+sound_speed;
            w_minus_c(it.index)=w-sound_speed;w_plus_c(it.index)=w+sound_speed;}}
    T dt_convect=max(u_minus_c.Max_Abs(),u_plus_c.Max_Abs())/grid.dX.x+max(v_minus_c.Max_Abs(),v_plus_c.Max_Abs())/grid.dX.y+max(w_minus_c.Max_Abs(),w_plus_c.Max_Abs())/grid.dX.z;
    return 1/dt_convect;
}              
//#####################################################################
template class EULER_3D<float>;
template class EULER_3D<double>;
