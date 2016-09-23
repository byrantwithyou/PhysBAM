//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Dynamics/Boundaries/BOUNDARY_REFLECTION_WATER.h>
namespace PhysBAM{
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_REFLECTION_WATER<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<T2,TV_INT>& u,ARRAY<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    assert(grid.Is_MAC_Grid() && phi && V);ARRAY<T2,TV_INT>::Put(u,u_ghost);
    RANGE<TV_INT> domain_indices=grid.Domain_Indices();
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int axis=0;axis<TV::m;axis++)for(int axis_side=0;axis_side<2;axis_side++){
        int side=2*axis+axis_side;
        if(use_extrapolation_mode && Constant_Extrapolation(side)) BOUNDARY<TV,T>::Fill_Single_Ghost_Region(grid,u_ghost,side,regions(side));
        else{ // either phi=phi_object for a wall, or no wall
            int inward_sign=axis_side==0?1:-1;T dx=grid.DX()[axis],half_dx=(T).5*dx;
            int cell_boundary=Boundary(side,regions(side)),face_boundary=cell_boundary+axis_side;
            int reflection_times_two=2*cell_boundary+(axis_side==0?-1:1);
            for(CELL_ITERATOR<TV> iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
                TV_INT boundary_cell=cell,boundary_face=cell;boundary_cell[axis]=cell_boundary;boundary_face[axis]=face_boundary;
                if((*phi)(boundary_cell) <= 0 && domain_indices.Lazy_Inside_Half_Open(boundary_cell) && inward_sign*V->Component(axis)(boundary_face) > tolerance){
                    TV_INT reflected_cell=cell;reflected_cell[axis]=reflection_times_two-cell[axis];
                    u_ghost(cell)=u_ghost(reflected_cell);}
                else u_ghost(cell)=u_ghost(boundary_cell);}}}
}  
}
