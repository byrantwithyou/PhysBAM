//#####################################################################
// Copyright 2009, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Function Use_Maccormack_Advection
//#####################################################################
// calculates the approximate area using Heaviside functions
template<class T_GRID> void LEVELSET_ADVECTION_UNIFORM<T_GRID>::
Use_Maccormack_Advection(const ARRAY<bool,TV_INT>& cell_mask)
{
    advection_maccormack=new ADVECTION_MACCORMACK_UNIFORM<T_GRID,T,ADVECTION<T_GRID,T> >(*advection,0,&cell_mask,0);
    Set_Custom_Advection(*advection_maccormack);
}
//#####################################################################
// Function Approximate_Negative_Material
//#####################################################################
// calculates the approximate area using Heaviside functions
template<class T_GRID> typename T_GRID::VECTOR_T::SCALAR LEVELSET_ADVECTION_UNIFORM<T_GRID>::
Approximate_Negative_Material(const T interface_thickness,const T time) const
{
    T_GRID& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    T_GRID node_grid=grid.Is_MAC_Grid()?grid.Get_Regular_Grid_At_MAC_Positions():grid;
    T interface_half_width=interface_thickness*grid.dX.Max()/2,volume=0;
    for(NODE_ITERATOR iterator(node_grid);iterator.Valid();iterator.Next()) volume+=LEVELSET_UTILITIES<T>::Heaviside(-phi(iterator.Node_Index()),interface_half_width);
    return volume*grid.Cell_Size();
}
//#####################################################################
// Function Approximate_Positive_Material
//#####################################################################
// calculates the approximate area using Heaviside functions
template<class T_GRID> typename T_GRID::VECTOR_T::SCALAR LEVELSET_ADVECTION_UNIFORM<T_GRID>::
Approximate_Positive_Material(const T interface_thickness,const T time) const
{
    T_GRID& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    T_GRID node_grid=grid.Is_MAC_Grid()?grid.Get_Regular_Grid_At_MAC_Positions():grid;
    T interface_half_width=interface_thickness*grid.dX.Max()/2,volume=0;
    for(NODE_ITERATOR iterator(node_grid);iterator.Valid();iterator.Next()) volume+=LEVELSET_UTILITIES<T>::Heaviside(phi(iterator.Node_Index()),interface_half_width);
    return volume*grid.Cell_Size();
}
//#####################################################################
// Function Euler_Step_Of_Reinitialization
//#####################################################################
template<class T_GRID> void LEVELSET_ADVECTION_UNIFORM<T_GRID>::
Euler_Step_Of_Reinitialization(const ARRAY<T,TV_INT>& sign_phi,const T dt,const T time)
{
    GRID<TV>& grid=levelset->grid;
    BOUNDARY_UNIFORM<GRID<TV>,T>* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    int ghost_cells=3;
    RANGE<TV_INT> domain(grid.Domain_Indices());
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));
    boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,ghost_cells);
    ARRAY<T,TV_INT> rhs(domain);

    for(int d=0;d<TV::m;d++){
        ARRAY<T,VECTOR<int,1> > phi_1d(-ghost_cells,grid.counts(d)+ghost_cells),phi_minus(0,grid.counts(d)),phi_plus(0,grid.counts(d));
        RANGE<VECTOR<int,TV::m-1> > domain_plane(VECTOR<int,TV::m-1>(),domain.max_corner.Remove_Index(d));
        for(RANGE_ITERATOR<TV::m-1> it(domain_plane);it.Valid();it.Next()){
            TV_INT i=it.index.Insert(0,d);
            for(i(d)=-ghost_cells;i(d)<grid.counts(d)+ghost_cells;i(d)++) phi_1d(i(d))=phi_ghost(i);
            if(reinitialization_spatial_order==5) HJ_WENO(grid.counts(d),grid.dX(d),phi_1d,phi_minus,phi_plus);
            else HJ_ENO(reinitialization_spatial_order,grid.counts(d),grid.dX(d),phi_1d,phi_minus,phi_plus);
            for(i(d)=0;i(d)<grid.counts(d);i(d)++)
                if(LEVELSET_UTILITIES<T>::Sign(phi(i))<0) rhs(i)=sqr(max(-phi_minus(i(d)),phi_plus(i(d)),(T)0));
                else rhs(i)=sqr(max(phi_minus(i(d)),-phi_plus(i(d)),(T)0));}}

    for(RANGE_ITERATOR<TV::m> it(domain);it.Valid();it.Next()){
        phi(it.index)-=dt*sign_phi(it.index)*(sqrt(rhs(it.index))-1);
        if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(it.index),phi(it.index))) phi(it.index)=LEVELSET_UTILITIES<T>::Sign(phi_ghost(it.index))*levelset->small_number*grid.min_dX;}

    boundary->Apply_Boundary_Condition(grid,phi,time); // time not incremented, pseudo-time
}
//#####################################################################
// Functions Reinitialize
//#####################################################################
template<class T_GRID> void LEVELSET_ADVECTION_UNIFORM<T_GRID>::
Reinitialize(const int time_steps,const T time)
{
    GRID<TV>& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;

    ARRAY<T,TV_INT> sign_phi(grid.Domain_Indices()); // smeared out sign function
    T epsilon=sqr(grid.dX.Max());
    TV_INT i;
    for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()) sign_phi(it.index)=phi(it.index)/sqrt(sqr(phi(it.index))+epsilon);

    T dt=reinitialization_cfl*grid.min_dX;
    RUNGEKUTTA<ARRAY<T,TV_INT> > rungekutta(phi);
    rungekutta.Set_Grid_And_Boundary_Condition(grid,*levelset->boundary);
    rungekutta.Set_Order(reinitialization_runge_kutta_order);
    rungekutta.Set_Time(time);
    rungekutta.Pseudo_Time();
    for(int k=0;k<time_steps;k++){
        rungekutta.Start(dt);
        for(int kk=0;kk<rungekutta.order;kk++){Euler_Step_Of_Reinitialization(sign_phi,dt,time);rungekutta.Main();}}
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void LEVELSET_ADVECTION_UNIFORM<T_GRID>::
Euler_Step(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocity,const T dt,const T time,const int number_of_ghost_cells)
{
    GRID<TV>& grid=levelset->grid;
    BOUNDARY_UNIFORM<GRID<TV>,T>* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    
    assert(grid.Is_MAC_Grid());
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(number_of_ghost_cells));
    boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);

    TV one_over_two_dX=(T).5*grid.one_over_dX,dphi;
    TV_INT off;
    off(TV::m-1)=1;
    for(int i=TV::m-1;i>0;i--) off(i-1)=off(i)*phi_ghost.counts(i);
    if(levelset->curvature_motion){ // do curvature first - based on phi^n
        bool curvature_defined=levelset->curvature!=0;levelset->Compute_Curvature(time);
        for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next()){
            int index=phi_ghost.Standard_Index(it.index);
            for(int i=0;i<TV::m;i++) dphi(i)=(phi_ghost.array(index+off(i))-phi_ghost.array(index-off(i)))*one_over_two_dX(i);
            phi(it.index)-=dt*levelset->sigma*(*levelset->curvature)(it.index)*dphi.Magnitude();}
        boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);
        if(!curvature_defined){delete levelset->curvature;levelset->curvature=0;}}

    advection->Update_Advection_Equation_Cell(grid,phi,phi_ghost,face_velocity,*boundary,dt,time);
    boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<float,1> > >;
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<float,2> > >;
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<double,1> > >;
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<double,2> > >;
template class LEVELSET_ADVECTION_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
