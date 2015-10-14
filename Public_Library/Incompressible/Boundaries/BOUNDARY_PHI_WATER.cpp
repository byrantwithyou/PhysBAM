//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_PHI_WATER
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BOUNDARY_PHI_WATER<TV>::
BOUNDARY_PHI_WATER(const TV_SIDES& constant_extrapolation)
    :sign(1),use_open_boundary_mode(false),open_boundary(6,true),V(0)
{
    Set_Constant_Extrapolation(constant_extrapolation);
    Set_Open_Boundary(false,false,false,false,false,false);
    Use_Extrapolation_Mode(false);
    Set_Tolerance();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BOUNDARY_PHI_WATER<TV>::
~BOUNDARY_PHI_WATER()
{}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV> void BOUNDARY_PHI_WATER<TV>::
Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_BASE& u,T_ARRAYS_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    assert(grid.Is_MAC_Grid() && V);T_ARRAYS_BASE::Put(u,u_ghost);
    RANGE<TV_INT> domain_indices=grid.Domain_Indices();
    VECTOR<RANGE<TV_INT>,2*TV::m> regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells); 
    for(int axis=0;axis<TV::m;axis++)for(int axis_side=0;axis_side<2;axis_side++){
        int side=2*axis+axis_side;
        Fill_Single_Ghost_Region_Threaded(regions(side), grid, u_ghost, side);}
}
//#####################################################################
// Function Fill_Single_Ghost_Region_Threaded
//#####################################################################
template<class TV> void BOUNDARY_PHI_WATER<TV>::
Fill_Single_Ghost_Region_Threaded(RANGE<TV_INT>& region,const GRID<TV>& grid,T_ARRAYS_BASE& u_ghost,const int side) const
{
    if(use_extrapolation_mode && Constant_Extrapolation(side)) BOUNDARY<TV,T>::Fill_Single_Ghost_Region(grid,u_ghost,side,region);
    else{ // either phi=phi_object for a wall, or no wall
        int axis_side=side%2;
        int axis=side/2;
        RANGE<TV_INT> domain_indices=grid.Domain_Indices();
        int inward_sign=axis_side==0?1:-1;T dx=grid.dX[axis],half_dx=(T).5*dx;
        int cell_boundary=Boundary(side,region),face_boundary=cell_boundary+axis_side;
        for(CELL_ITERATOR<TV> iterator(grid,region);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            TV_INT boundary_cell=cell,boundary_face=cell;boundary_cell[axis]=cell_boundary;boundary_face[axis]=face_boundary;
            if(use_open_boundary_mode&&open_boundary(side)) u_ghost(cell)=callbacks->Open_Water_Phi(iterator.Location(),0);
            else if(sign*u_ghost(boundary_cell) <= 0 && domain_indices.Lazy_Inside_Half_Open(boundary_cell) && inward_sign*V->Component(axis)(boundary_face) > tolerance){
                T distance_to_boundary=inward_sign*dx*(cell_boundary-cell[axis]);
                u_ghost(cell)=sign*(distance_to_boundary+max(-half_dx,sign*u_ghost(boundary_cell)));} // distance to wall
            else u_ghost(cell)=u_ghost(boundary_cell);}}
}
//#####################################################################
template class BOUNDARY_PHI_WATER<VECTOR<float,1> >;
template class BOUNDARY_PHI_WATER<VECTOR<float,2> >;
template class BOUNDARY_PHI_WATER<VECTOR<float,3> >;
template class BOUNDARY_PHI_WATER<VECTOR<double,1> >;
template class BOUNDARY_PHI_WATER<VECTOR<double,2> >;
template class BOUNDARY_PHI_WATER<VECTOR<double,3> >;
}
