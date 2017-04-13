//#####################################################################
// Copyright 2002-2004, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER  
//##################################################################### 
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_EIGENSYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SHALLOW_WATER<TV>::
SHALLOW_WATER(GRID<TV>& grid,ARRAY<TV_DIMENSION,TV_INT>& U,
    T gravity,T min_height)
    :gravity(gravity),min_height(min_height),grid(grid),U(U)
{
    boundary=&boundary_default;
    conservation=&conservation_default;
    for(int i=0;i<TV::m;i++)
        eigensystems(i)=new SHALLOW_WATER_EIGENSYSTEM<TV>(gravity,i);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SHALLOW_WATER<TV>::
~SHALLOW_WATER()
{
    for(int i=0;i<TV::m;i++)
        delete eigensystems(i);
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class TV> void SHALLOW_WATER<TV>::
Euler_Step(const T dt,const T time)
{  
    int ghost_cells=3;

    // make sure things'll work in conservation law solver
    TV_DIMENSION min_u;
    min_u(0)=min_height;
    for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        TV_DIMENSION& u=U(it.index);
        if(u(0)<min_height) u=min_u;}

    ARRAY<TV_DIMENSION,TV_INT> U_ghost(grid.Node_Indices(ghost_cells));
    boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,ghost_cells);

    ARRAY<bool,TV_INT> psi(grid.Node_Indices(ghost_cells),true,true); // no cut out grids
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N(grid.Get_MAC_Grid_At_Regular_Positions());
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities(grid.Get_MAC_Grid_At_Regular_Positions());
    conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystems,eigensystems,psi_N,face_velocities);
    boundary->Apply_Boundary_Condition(grid,U,time+dt); 
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR SHALLOW_WATER<TV>::
CFL()
{
    TV max_speed;
    for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        TV_DIMENSION UU=U(it.index);
        T h=UU(0),celerity=sqrt(gravity*h);
        TV u=UU.template Slice<1,TV::m>()/h;
        max_speed=max(max_speed,abs(u)+celerity);}
    T dt_convect=max_speed.Dot(grid.one_over_dX);
    return 1/dt_convect;
}
//#####################################################################
template class SHALLOW_WATER<VECTOR<float,1> >;
template class SHALLOW_WATER<VECTOR<float,2> >;
template class SHALLOW_WATER<VECTOR<double,1> >;
template class SHALLOW_WATER<VECTOR<double,2> >;
}
