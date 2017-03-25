//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Jon Gretarsson, Eran Guendelman, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Poisson/POISSON_UNIFORM.h>
using namespace PhysBAM;
template<class TV> POISSON_UNIFORM<TV>::
POISSON_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input)
    :POISSON<T>(multiphase_input),LAPLACE_UNIFORM<TV>(grid_input,u_input,false,enforce_compatibility_input)
{
    Initialize_Grid(grid_input);
}
template<class TV> POISSON_UNIFORM<TV>::
~POISSON_UNIFORM()
{
}
//#####################################################################
// Function Find_Variable_beta
//#####################################################################
// only set up for Dirichlet boundary conditions - doesn't work for jump conditons yet 
template<class TV> void POISSON_UNIFORM<TV>::
Find_Variable_beta()
{
    if(beta_given_on_faces) return; // beta_right, beta_top, beta_back already set
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) 
        beta_face.Component(iterator.Axis())(iterator.Face_Index())=(T).5*(variable_beta(iterator.First_Cell_Index())+variable_beta(iterator.Second_Cell_Index())); 
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class TV> void POISSON_UNIFORM<TV>::
Find_A_Part_Two(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    TV one_over_dx2=Inverse(grid.dX*grid.dX);
    TV_INT grid_counts=grid.counts;
    if(use_weighted_divergence)
        for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){
            int color=filled_region_colors(iterator.Cell_Index());
            if(color!=-2 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){const TV_INT& cell_index=iterator.Cell_Index();
                int matrix_index=cell_index_to_matrix_index(cell_index);
                SPARSE_MATRIX_FLAT_MXN<T>& A=A_array(filled_region_colors(cell_index));ARRAY<T>& b=b_array(filled_region_colors(cell_index));b(matrix_index)=f(cell_index);
                T diagonal=0;
                for(int axis=0;axis<TV::m;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
                    if(filled_region_colors.Valid_Index(cell_index-offset)){
                        if(!psi_N.Component(axis)(cell_index)){
                            T element=divergence_face_weights.Component(axis)(cell_index)*beta_face.Component(axis)(cell_index)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside_Half_Open(cell_index-offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index-offset;
                                int axis_periodic_cell=wrap(periodic_offset_cell[axis],grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index-offset)) b(matrix_index)-=element*u(cell_index-offset);
                            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),element);}}
                    if(filled_region_colors.Valid_Index(cell_index+offset)){
                        if(!psi_N.Component(axis)(cell_index+offset)){
                            T element=divergence_face_weights.Component(axis)(cell_index+offset)*beta_face.Component(axis)(cell_index+offset)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside_Half_Open(cell_index+offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index+offset;
                                int axis_periodic_cell=wrap(periodic_offset_cell[axis],grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index+offset)) b(matrix_index)-=element*u(cell_index+offset);
                            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),element);}}}
                A.Set_Element(matrix_index,matrix_index,diagonal);}} // set diagonal and right hand side
    else
        for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){
            int color=filled_region_colors(iterator.Cell_Index());
            if(color!=-2 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){const TV_INT& cell_index=iterator.Cell_Index();
                int matrix_index=cell_index_to_matrix_index(cell_index);
                SPARSE_MATRIX_FLAT_MXN<T>& A=A_array(filled_region_colors(cell_index));ARRAY<T>& b=b_array(filled_region_colors(cell_index));b(matrix_index)=f(cell_index);
                T diagonal=0;
                for(int axis=0;axis<TV::m;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
                    if(filled_region_colors.Valid_Index(cell_index-offset)){
                        if(!psi_N.Component(axis)(cell_index)){T element=beta_face.Component(axis)(cell_index)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside_Half_Open(cell_index-offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index-offset;
                                int axis_periodic_cell=wrap(periodic_offset_cell[axis],grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index-offset)) b(matrix_index)-=element*u(cell_index-offset);
                            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),element);}}
                    if(filled_region_colors.Valid_Index(cell_index+offset)){
                        if(!psi_N.Component(axis)(cell_index+offset)){T element=beta_face.Component(axis)(cell_index+offset)*one_over_dx2[axis];
                            diagonal-=element;
                            if(grid.Domain_Indices().Lazy_Outside_Half_Open(cell_index+offset) && periodic_boundary[axis]){
                                TV_INT periodic_offset_cell=cell_index+offset;
                                int axis_periodic_cell=wrap(periodic_offset_cell[axis],grid_counts[axis]);
                                periodic_offset_cell[axis]=axis_periodic_cell;
                                A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),element);}
                            else if(psi_D(cell_index+offset)) b(matrix_index)-=element*u(cell_index+offset);
                            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),element);}}}
                A.Set_Element(matrix_index,matrix_index,diagonal);}} // set diagonal and right hand side
}
template<class TV> void POISSON_UNIFORM<TV>::
Initialize_Grid(const GRID<TV>& grid_input)
{
    LAPLACE_UNIFORM<TV>::Initialize_Grid(grid_input); // laplace will call Set_Up_Second_Order_Cut_Cell_Method to resize those arrays
    beta_face.Resize(grid,1); // TODO: check if can get rid of extra ghost layer
    if(use_variable_beta)variable_beta.Resize(grid.Domain_Indices(1));
}
//#####################################################################
namespace PhysBAM{
template class POISSON_UNIFORM<VECTOR<float,1> >;
template class POISSON_UNIFORM<VECTOR<float,2> >;
template class POISSON_UNIFORM<VECTOR<float,3> >;
template class POISSON_UNIFORM<VECTOR<double,1> >;
template class POISSON_UNIFORM<VECTOR<double,2> >;
template class POISSON_UNIFORM<VECTOR<double,3> >;
}
