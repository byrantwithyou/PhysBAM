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
