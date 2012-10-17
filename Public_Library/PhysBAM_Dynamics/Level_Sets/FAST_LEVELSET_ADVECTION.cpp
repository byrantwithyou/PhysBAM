//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/FAST_LEVELSET_ADVECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Local_WENO_Advect
//#####################################################################
// phi is (-2,m+3), u, distance and u_phix are [0,m)
template<class T> static void
Local_WENO_Advect(const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& u,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& u_phix,const T half_band_width)
{
    T epsilon=(T)1e-6*sqr(dx); // 1e-6 works since phi is a distance function - sqr(dx) since undivided differences are used
    T one_over_dx=1/dx;
    for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){ // one_over_dx since undivided differences are used
        if(u(i) > 0) u_phix(i)=u(i)*ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i-2)-phi(i-3),phi(i-1)-phi(i-2),phi(i)-phi(i-1),phi(i+1)-phi(i),phi(i+2)-phi(i+1),epsilon)*one_over_dx;
        else u_phix(i)=u(i)*ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i+3)-phi(i+2),phi(i+2)-phi(i+1),phi(i+1)-phi(i),phi(i)-phi(i-1),phi(i-1)-phi(i-2),epsilon)*one_over_dx;}
}
//#####################################################################
// Function Local_ENO_Advect
//#####################################################################
// order = 1, 2 or 3, phi is (-2,m_3), u, distance and u_phix are [0,m)
template<class T> static void
Local_ENO_Advect(const int order,const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& u,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& u_phix,const T half_band_width)
{
    T one_over_dx=1/dx;
    if(order == 1){for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){
        if(u(i) > 0) u_phix(i)=u(i)*(phi(i)-phi(i-1))*one_over_dx;
        else u_phix(i)=u(i)*(phi(i+1)-phi(i))*one_over_dx;}}
    else if(order == 2){for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){
        if(u(i) > 0) u_phix(i)=u(i)*(phi(i)-phi(i-1)+(T).5*minmag(phi(i)-2*phi(i-1)+phi(i-2),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;
        else u_phix(i)=u(i)*(phi(i+1)-phi(i)-(T).5*minmag(phi(i+2)-2*phi(i+1)+phi(i),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;}}
    else if(order == 3){for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){
        if(u(i) > 0){
            T D2_left=phi(i)-2*phi(i-1)+phi(i-2),D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) 2*dx*dx
            if(abs(D2_left) <= abs(D2_right))
                u_phix(i)=u(i)*(phi(i)-phi(i-1)+(T).5*D2_left+(T)one_third*minmag(phi(i)-3*(phi(i-1)-phi(i-2))-phi(i-3),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;
            else u_phix(i)=u(i)*(phi(i)-phi(i-1)+(T).5*D2_right-(T)one_sixth*minmag(phi(i+2)-3*(phi(i+1)-phi(i))-phi(i-1),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;}
        else{
            T D2_left=phi(i+2)-2*phi(i+1)+phi(i),D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) -2*dx*dx
            if(abs(D2_left) <= abs(D2_right))
                u_phix(i)=u(i)*(phi(i+1)-phi(i)-(T).5*D2_left+(T)one_third*minmag(phi(i+3)-3*(phi(i+2)-phi(i+1))+phi(i),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;
            else u_phix(i)=u(i)*(phi(i+1)-phi(i)-(T).5*D2_right-(T)one_sixth*minmag(phi(i+1)-3*(phi(i)-phi(i-1))+phi(i-2),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;}}}
}
//#####################################################################
// Function Euler_Step_High_Order_Helper
//#####################################################################
template<class T> static void
Euler_Step_High_Order_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAY<VECTOR<T,1> ,VECTOR<int,1> >& V,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& phi_ghost,ARRAY<T,VECTOR<int,1> >& rhs,const int spatial_order,
    const T half_band_width)
{
    int m=grid.counts.x;T dx=grid.dX.x;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(-ghost_cells,m+ghost_cells),u_1d(0,m),distance_1d_x(0,m);
    for(int i=-ghost_cells;i<m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i);
    for(int i=0;i<m;i++){u_1d(i)=V(i).x;distance_1d_x(i)=phi(i);}
    if(spatial_order == 5) Local_WENO_Advect(m,dx,phi_1d_x,u_1d,distance_1d_x,rhs,half_band_width);
    else Local_ENO_Advect(spatial_order,m,dx,phi_1d_x,u_1d,distance_1d_x,rhs,half_band_width);
}
template<class T> static void
Euler_Step_High_Order_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V,const ARRAY<T,VECTOR<int,2> >& phi,const ARRAY<T,VECTOR<int,2> >& phi_ghost,ARRAY<T,VECTOR<int,2> >& rhs,const int spatial_order,
    const T half_band_width)
{
    int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(-ghost_cells,m+ghost_cells),u_1d(0,m),distance_1d_x(0,m),u_phix_1d(0,m);
    for(int j=0;j<n;j++){
        for(int i=-ghost_cells;i<m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i,j);
        for(int i=0;i<m;i++){u_1d(i)=V(i,j).x;distance_1d_x(i)=phi(i,j);}
        if(spatial_order == 5) Local_WENO_Advect(m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        for(int i=0;i<m;i++) rhs(i,j)=u_phix_1d(i);}
    ARRAY<T,VECTOR<int,1> > phi_1d_y(-ghost_cells,n+ghost_cells),v_1d(0,n),distance_1d_y(0,n),v_phiy_1d(0,n);
    for(int i=0;i<m;i++){
        for(int j=-ghost_cells;j<n+ghost_cells;j++) phi_1d_y(j)=phi_ghost(i,j);
        for(int j=0;j<n;j++){v_1d(j)=V(i,j).y;distance_1d_y(j)=phi(i,j);}
        if(spatial_order == 5) Local_WENO_Advect(n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        for(int j=0;j<n;j++) rhs(i,j)+=v_phiy_1d(j);}
}
template<class T> static void
Euler_Step_High_Order_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& V,const ARRAY<T,VECTOR<int,3> >& phi,const ARRAY<T,VECTOR<int,3> >& phi_ghost,ARRAY<T,VECTOR<int,3> >& rhs,const int spatial_order,
    const T half_band_width)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;T dx=grid.dX.x,dy=grid.dX.y,dz=grid.dX.z;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(-ghost_cells,m+ghost_cells),u_1d(0,m),distance_1d_x(0,m),u_phix_1d(0,m);
    for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++){
        for(int i=-ghost_cells;i<m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i,j,ij);
        for(int i=0;i<m;i++){u_1d(i)=V(i,j,ij).x;distance_1d_x(i)=phi(i,j,ij);}
        if(spatial_order == 5) Local_WENO_Advect(m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        for(int i=0;i<m;i++) rhs(i,j,ij)=u_phix_1d(i);}
    ARRAY<T,VECTOR<int,1> > phi_1d_y(-ghost_cells,n+ghost_cells),v_1d(0,n),distance_1d_y(0,n),v_phiy_1d(0,n);
    for(int i=0;i<m;i++) for(int ij=0;ij<mn;ij++){
        for(int j=-ghost_cells;j<n+ghost_cells;j++) phi_1d_y(j)=phi_ghost(i,j,ij);
        for(int j=0;j<n;j++){v_1d(j)=V(i,j,ij).y;distance_1d_y(j)=phi(i,j,ij);}
        if(spatial_order == 5) Local_WENO_Advect(n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        for(int j=0;j<n;j++) rhs(i,j,ij)+=v_phiy_1d(j);}
    ARRAY<T,VECTOR<int,1> > phi_1d_z(-ghost_cells,mn+ghost_cells),w_1d(0,mn),distance_1d_z(0,mn),w_phiz_1d(0,mn);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        for(int ij=-ghost_cells;ij<mn+ghost_cells;ij++) phi_1d_z(ij)=phi_ghost(i,j,ij);
        for(int ij=0;ij<mn;ij++){w_1d(ij)=V(i,j,ij).z;distance_1d_z(ij)=phi(i,j,ij);}
        if(spatial_order == 5) Local_WENO_Advect(mn,dz,phi_1d_z,w_1d,distance_1d_z,w_phiz_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,mn,dz,phi_1d_z,w_1d,distance_1d_z,w_phiz_1d,half_band_width);
        for(int ij=0;ij<mn;ij++) rhs(i,j,ij)+=w_phiz_1d(ij);}
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void FAST_LEVELSET_ADVECTION<T_GRID>::
Euler_Step(const ARRAY<TV,TV_INT>& V,const T dt,const T time,const int number_of_ghost_cells)
{
    FAST_LEVELSET<GRID<TV> >* fl=(FAST_LEVELSET<GRID<TV> >*)levelset;
    T_GRID& grid=fl->grid;
    T_ARRAYS_SCALAR& phi=fl->phi;
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(number_of_ghost_cells));fl->boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);

    if(local_semi_lagrangian_advection){
        LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(phi_ghost(cell)) <= fl->half_band_width)
            phi(cell)=interpolation.Clamped_To_Array(grid,phi_ghost,iterator.Location()-dt*V(cell));}}
    else if(local_advection_spatial_order){
        T_ARRAYS_SCALAR rhs(grid.Domain_Indices());
        Euler_Step_High_Order_Helper(grid,V,phi,phi_ghost,rhs,local_advection_spatial_order,fl->half_band_width);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(phi(cell)) <= fl->half_band_width) phi(cell)-=dt*rhs(cell);}}
    else // use the advection routine in the level set base class
        advection->Update_Advection_Equation_Node(grid,phi,phi_ghost,V,*fl->boundary,dt,time);

    fl->boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void FAST_LEVELSET_ADVECTION<T_GRID>::
Euler_Step(const T_FACE_ARRAYS_SCALAR& V,const T dt,const T time,const int number_of_ghost_cells)
{
    FAST_LEVELSET<GRID<TV> >* fl=(FAST_LEVELSET<GRID<TV> >*)levelset;
    T_GRID& grid=fl->grid;
    T_ARRAYS_SCALAR& phi=fl->phi;
    
    assert(grid.Is_MAC_Grid() && advection); // for now use advection in base class
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(number_of_ghost_cells));fl->boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);
    advection->Update_Advection_Equation_Cell(grid,phi,phi_ghost,V,*fl->boundary,dt,time);
    fl->boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}
//#####################################################################
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<float,1> > >;
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<float,2> > >;
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<double,1> > >;
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<double,2> > >;
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<double,3> > >;
#endif
