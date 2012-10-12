//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Local_ENO_Reinitialize
//#####################################################################
// order = 1, 2 or 3, phi is (-2,m_3), phix_minus and phix_plus are [0,m)
template<class T> static void
Local_ENO_Reinitialize(const int order,const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus,const T half_band_width)
{
    T one_over_dx=1/dx;
    if(order == 1){for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){
        phix_minus(i)=(phi(i)-phi(i-1))*one_over_dx;
        phix_plus(i)=(phi(i+1)-phi(i))*one_over_dx;}}
    else if(order == 2){for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){
        phix_minus(i)=(phi(i)-phi(i-1)+(T).5*minmag(phi(i)-2*phi(i-1)+phi(i-2),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;
        phix_plus(i)=(phi(i+1)-phi(i)-(T).5*minmag(phi(i+2)-2*phi(i+1)+phi(i),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;}}
    else if(order == 3){for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){
        T D2_left=phi(i)-2*phi(i-1)+phi(i-2),D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) 2*dx*dx
        if(abs(D2_left) <= abs(D2_right))
            phix_minus(i)=(phi(i)-phi(i-1)+(T).5*D2_left+(T)one_third*minmag(phi(i)-3*(phi(i-1)-phi(i-2))-phi(i-3),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;
        else phix_minus(i)=(phi(i)-phi(i-1)+(T).5*D2_right-(T)one_sixth*minmag(phi(i+2)-3*(phi(i+1)-phi(i))-phi(i-1),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;
        D2_left=phi(i+2)-2*phi(i+1)+phi(i);D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) -2*dx*dx
        if(abs(D2_left) <= abs(D2_right))
            phix_plus(i)=(phi(i+1)-phi(i)-(T).5*D2_left+(T)one_third*minmag(phi(i+3)-3*(phi(i+2)-phi(i+1))+phi(i),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;
        else phix_plus(i)=(phi(i+1)-phi(i)-(T).5*D2_right-(T)one_sixth*minmag(phi(i+1)-3*(phi(i)-phi(i-1))+phi(i-2),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;}}
}
//#####################################################################
// Function Local_WENO_Reinitialize
//#####################################################################
// phi is (-2,m+3), distance and phix_minus and phix_plus are [0,m)
template<class T> static void
Local_WENO_Reinitialize(const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus,const T half_band_width)
{
    T epsilon=(T)1e-6*sqr(dx); // 1e-6 works since phi is a distance function - sqr(dx) since undivided differences are used
    T one_over_dx=1/dx;
    for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){ // one_over_dx since undivided differences are used
        phix_minus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i-2)-phi(i-3),phi(i-1)-phi(i-2),phi(i)-phi(i-1),phi(i+1)-phi(i),phi(i+2)-phi(i+1),epsilon)*one_over_dx;
        phix_plus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i+3)-phi(i+2),phi(i+2)-phi(i+1),phi(i+1)-phi(i),phi(i)-phi(i-1),phi(i-1)-phi(i-2),epsilon)*one_over_dx;}
}
//#####################################################################
// Functions Euler_Step_Of_Reinitialization
//#####################################################################
template<class T,class TV,class TV_INT> static void
Euler_Step_Of_Reinitialization(FAST_LEVELSET<GRID<TV> >* levelset,const ARRAY<T,TV_INT>& signed_distance,const ARRAY<T,TV_INT>& sign_phi,T dt,T time,T half_band_width,int spatial_order)
{
    GRID<TV>& grid=levelset->grid;
    ARRAY<T,TV_INT>& phi=levelset->phi;
    
    int ghost_cells=3;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));
    levelset->boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,ghost_cells);
    ARRAY<T,TV_INT> rhs(grid.Domain_Indices());

    for(int d=0;d<TV::m;d++){
        int m=grid.counts(d);
        ARRAY<T,VECTOR<int,1> > phi_1d(-ghost_cells,m+3),distance_1d(0,m),phi_minus(0,m),phi_plus(0,m);
        RANGE<TV_INT> range(TV_INT(),grid.counts);
        range.max_corner(d)=1;
        for(RANGE_ITERATOR<TV::m> it(range);it.Valid();it.Next()){
            TV_INT index(it.index);
            int& i=index(d);
            for(i=-ghost_cells;i<m+3;i++) phi_1d(i)=phi_ghost(index);
            for(i=0;i<m;i++) distance_1d(i)=signed_distance(index);
            if(spatial_order==5)
                Local_WENO_Reinitialize(m,grid.dX(d),phi_1d,distance_1d,phi_minus,phi_plus,half_band_width);
            else Local_ENO_Reinitialize(spatial_order,m,grid.dX(d),phi_1d,distance_1d,phi_minus,phi_plus,half_band_width);
            for(i=0;i<m;i++)
                if(abs(signed_distance(index))<=half_band_width){
                    if(LEVELSET_UTILITIES<T>::Sign(phi(index)) < 0) rhs(index)+=sqr(max(-phi_minus(i),phi_plus(i),(T)0));
                    else rhs(index)+=sqr(max(phi_minus(i),-phi_plus(i),(T)0));}}}

    T min_DX=grid.dX.Min();
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(signed_distance(cell)) <= half_band_width){
        phi(cell)-=dt*sign_phi(cell)*(sqrt(rhs(cell))-1);
        if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(cell),phi(cell))) phi(cell)=LEVELSET_UTILITIES<T>::Sign(phi_ghost(cell))*levelset->small_number*min_DX;}}

    levelset->boundary->Apply_Boundary_Condition(grid,phi,time); // time not incremented - pseudo-time
}
//#####################################################################
// Functions Reinitialize
//#####################################################################
// extra_band=grid.dX.Max()*(1+min(3,local_advection_spatial_order))
namespace PhysBAM{
template<class T,class TV> void
Reinitialize(FAST_LEVELSET<GRID<TV> >* levelset,int time_steps,T time,T half_band_width,T extra_band,T cfl,int temporal_order,int spatial_order)
{
    typedef VECTOR<int,TV::m> TV_INT;
    GRID<TV>& grid=levelset->grid;
    ARRAY<T,TV_INT>& phi=levelset->phi;
    
    T large_band=half_band_width+extra_band;
    ARRAY<T,TV_INT> signed_distance(grid.Domain_Indices());
    levelset->Get_Signed_Distance_Using_FMM(signed_distance,time,large_band);

    ARRAY<T,TV_INT> sign_phi(grid.Domain_Indices()); // smeared out sign function
    T epsilon=sqr(grid.dX.Max());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        sign_phi(cell)=phi(cell)/sqrt(sqr(phi(cell))+epsilon);}

    T dt=cfl*grid.dX.Min();
    RUNGEKUTTA<ARRAY<T,TV_INT> > rungekutta(phi);
    rungekutta.Set_Grid_And_Boundary_Condition(grid,*levelset->boundary);
    rungekutta.Set_Order(temporal_order);
    rungekutta.Set_Time(time);
    rungekutta.Pseudo_Time();
    for(int k=0;k<time_steps;k++){
        rungekutta.Start(dt);
        for(int kk=0;kk<rungekutta.order;kk++){
            Euler_Step_Of_Reinitialization(levelset,signed_distance,sign_phi,dt,time,half_band_width,spatial_order);
            rungekutta.Main();}
        }

    T min_DX=grid.dX.Min();
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(signed_distance(cell)) <= large_band){
        if(abs(signed_distance(cell)) > half_band_width) phi(cell)=signed_distance(cell); // outer band - use the FMM solution
        else if(abs(signed_distance(cell)-phi(cell)) > min_DX) phi(cell)=signed_distance(cell);}} // inner band - use FMM if errors look big
}
template void Reinitialize<float,VECTOR<float,1> >(FAST_LEVELSET<GRID<VECTOR<float,1> > >*,int,float,float,float,float,int,int);
template void Reinitialize<float,VECTOR<float,2> >(FAST_LEVELSET<GRID<VECTOR<float,2> > >*,int,float,float,float,float,int,int);
template void Reinitialize<float,VECTOR<float,3> >(FAST_LEVELSET<GRID<VECTOR<float,3> > >*,int,float,float,float,float,int,int);
template void Reinitialize<double,VECTOR<double,1> >(FAST_LEVELSET<GRID<VECTOR<double,1> > >*,int,double,double,double,double,int,int);
template void Reinitialize<double,VECTOR<double,2> >(FAST_LEVELSET<GRID<VECTOR<double,2> > >*,int,double,double,double,double,int,int);
template void Reinitialize<double,VECTOR<double,3> >(FAST_LEVELSET<GRID<VECTOR<double,3> > >*,int,double,double,double,double,int,int);
}
