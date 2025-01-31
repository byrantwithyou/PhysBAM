//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_COLLIDABLE_UNIFORM  
//#####################################################################
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/LAPLACE_COLLIDABLE_UNIFORM.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LAPLACE_COLLIDABLE_UNIFORM<TV>::
LAPLACE_COLLIDABLE_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input)
    :BASE(grid_input,u_input,initialize_grid,enforce_compatibility_input),levelset_default(new LEVELSET<TV>(grid,phi_default))
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LAPLACE_COLLIDABLE_UNIFORM<TV>::
LAPLACE_COLLIDABLE_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,LEVELSET<TV>& cell_centered_levelset,const bool initialize_grid,const bool multiphase_input,
    const bool enforce_compatibility_input)
    :BASE(grid_input,u_input,initialize_grid,enforce_compatibility_input),levelset_default(new LEVELSET<TV>(grid,phi_default))
{
    levelset=&cell_centered_levelset;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LAPLACE_COLLIDABLE_UNIFORM<TV>::
~LAPLACE_COLLIDABLE_UNIFORM()
{
    delete levelset_default;
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class TV> void LAPLACE_COLLIDABLE_UNIFORM<TV>::
Find_A(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    BASE::Find_A(A_array,b_array,filled_region_cell_count,cell_index_to_matrix_index);
    if(second_order_cut_cell_method) Apply_Second_Order_Cut_Cell_Method(A_array,b_array,cell_index_to_matrix_index);
}
//#####################################################################
// Function Apply_Second_Order_Cut_Cell_Method
//#####################################################################
template<class TV> void LAPLACE_COLLIDABLE_UNIFORM<TV>::
Apply_Second_Order_Cut_Cell_Method(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    assert(levelset);
    TV minus_one_over_dx_squared=(T)-1*Inverse(grid.dX*grid.dX);
    for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::INTERIOR_REGION);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index(),cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
        if(!psi_N.Component(axis)(face) && LEVELSET_UTILITIES<T>::Interface(levelset->phi(cell1),levelset->phi(cell2))){
            T theta=LEVELSET_UTILITIES<T>::Theta(levelset->phi(cell1),levelset->phi(cell2));
            if(psi_D(cell1) && !psi_D(cell2) && grid.Domain_Indices(1).Lazy_Inside_Half_Open(cell2)){ // interface is to the negative side of second cell
                int color=filled_region_colors(cell2);int matrix_index=cell_index_to_matrix_index(cell2);
                T A_right_i=minus_one_over_dx_squared[axis];A_array(color).Add_Element(matrix_index,matrix_index,-A_right_i);b_array(color)(matrix_index)-=A_right_i*u(cell1);
                A_right_i/=max((1-theta),second_order_cut_cell_threshold);
                A_array(color).Add_Element(matrix_index,matrix_index,A_right_i);b_array(color)(matrix_index)+=A_right_i*u_interface.Component(axis)(face);}  
            else if(psi_D(cell2) && !psi_D(cell1) && grid.Domain_Indices(1).Lazy_Inside_Half_Open(cell1)){ // interface is to the positive side of first cell
                int color=filled_region_colors(cell1);int matrix_index=cell_index_to_matrix_index(cell1);
                T A_right_i=minus_one_over_dx_squared[axis];A_array(color).Add_Element(matrix_index,matrix_index,-A_right_i);b_array(color)(matrix_index)-=A_right_i*u(cell2);
                A_right_i/=max(theta,second_order_cut_cell_threshold);
                A_array(color).Add_Element(matrix_index,matrix_index,A_right_i);b_array(color)(matrix_index)+=A_right_i*u_interface.Component(axis)(face);}}}
}
template<class TV> void LAPLACE_COLLIDABLE_UNIFORM<TV>::
Initialize_Grid(const GRID<TV>& mac_grid_input)
{
    assert(mac_grid_input.DX()==TV() || mac_grid_input.Is_MAC_Grid());
    BASE::Initialize_Grid(mac_grid_input);
    if(levelset == levelset_default) phi_default.Resize(grid.Domain_Indices(1));
    Set_Up_Second_Order_Cut_Cell_Method(second_order_cut_cell_method);
}
template<class TV> void LAPLACE_COLLIDABLE_UNIFORM<TV>::
Set_Up_Second_Order_Cut_Cell_Method(const bool use_second_order_cut_cell_method_input)
{
    second_order_cut_cell_method=use_second_order_cut_cell_method_input;
    if(second_order_cut_cell_method) u_interface.Resize(grid);
    else u_interface.Clean_Memory();
}
//#####################################################################
template class LAPLACE_COLLIDABLE_UNIFORM<VECTOR<float,1> >;
template class LAPLACE_COLLIDABLE_UNIFORM<VECTOR<float,2> >;
template class LAPLACE_COLLIDABLE_UNIFORM<VECTOR<float,3> >;
template class LAPLACE_COLLIDABLE_UNIFORM<VECTOR<double,1> >;
template class LAPLACE_COLLIDABLE_UNIFORM<VECTOR<double,2> >;
template class LAPLACE_COLLIDABLE_UNIFORM<VECTOR<double,3> >;
}
