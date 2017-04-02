//#####################################################################
// Copyright 2002-2004, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D  
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T> void SHALLOW_WATER_2D<T>::
Euler_Step(const T dt,const T time)
{  
    int m=grid.counts.x,n=grid.counts.y;
    int ghost_cells=3;

    // make sure things'll work in conservation law solver
    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++) if(U(i,j)(0) < min_height){U(i,j)(0)=min_height;U(i,j)(1)=U(i,j)(2)=0;}

    ARRAY<TV_DIMENSION,VECTOR<int,2> > U_ghost(-ghost_cells,m+ghost_cells,-ghost_cells,n+ghost_cells);boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);

    ARRAY<bool,VECTOR<int,2> > psi(0,m,0,n);psi.Fill(true); // no cut out grids
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    VECTOR<EIGENSYSTEM<T,3>*,2> eigensystem(&eigensystem_F,&eigensystem_G);
    conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystem,eigensystem,psi_N,face_velocities);

    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T> T SHALLOW_WATER_2D<T>::
CFL()
{
    T max_x_speed=0,max_y_speed=0;
    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++){
        T one_over_h=1/U(i,j)(0);
        T u=U(i,j)(1)*one_over_h,v=U(i,j)(2)*one_over_h,celerity=sqrt(gravity*U(i,j)(0));
        max_x_speed=max(max_x_speed,abs(u)+celerity);
        max_y_speed=max(max_y_speed,abs(v)+celerity);}
    T dt_convect=max_x_speed*grid.one_over_dX.x+max_y_speed*grid.one_over_dX.y;
    return 1/dt_convect;
}
//#####################################################################
template class SHALLOW_WATER_2D<float>;
template class SHALLOW_WATER_2D<double>;
