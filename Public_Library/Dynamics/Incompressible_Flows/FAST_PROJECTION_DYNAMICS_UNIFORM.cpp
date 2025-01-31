//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Dynamics/Incompressible_Flows/FAST_PROJECTION_DYNAMICS_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FAST_PROJECTION_DYNAMICS_UNIFORM<TV>::
FAST_PROJECTION_DYNAMICS_UNIFORM(const int scale,const bool flame_input,const bool multiphase,const bool use_variable_beta,const bool use_poisson)
    :PROJECTION_DYNAMICS_UNIFORM<TV>(GRID<TV>(TV_INT::All_Ones_Vector()*scale,RANGE<TV>(TV(),TV::All_Ones_Vector()),true),flame_input,multiphase,use_variable_beta,use_poisson),
    A(*new SPARSE_MATRIX_FLAT_MXN<T>)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FAST_PROJECTION_DYNAMICS_UNIFORM<TV>::
FAST_PROJECTION_DYNAMICS_UNIFORM(const int scale,LEVELSET<TV>& levelset_input)
    :PROJECTION_DYNAMICS_UNIFORM<TV>(GRID<TV>(TV_INT::All_Ones_Vector()*scale,RANGE<TV>(TV(),TV::All_Ones_Vector()),true),levelset_input),
    A(*new SPARSE_MATRIX_FLAT_MXN<T>)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FAST_PROJECTION_DYNAMICS_UNIFORM<TV>::
~FAST_PROJECTION_DYNAMICS_UNIFORM()
{
    delete &A;
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class TV> void FAST_PROJECTION_DYNAMICS_UNIFORM<TV>::
Initialize_Grid(const GRID<TV>& mac_grid)
{
    BASE::Initialize_Grid(mac_grid);
    int scale=mac_grid.counts.x;
    elliptic_solver->Set_Neumann_Outer_Boundaries();
    int number_of_elements=TV::m==2?scale*scale:scale*scale*scale;
    cell_index_to_matrix_index.Resize(mac_grid.Domain_Indices());
    cell_index_to_matrix_index.Fill(-1);
    matrix_index_to_cell_index.Resize(number_of_elements);
    b.Resize(number_of_elements);
    ARRAY<int> row_counts;
    row_counts.Resize(number_of_elements);
    int count=0;
    for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        int matrix_index;
        cell_index_to_matrix_index(cell_index)=matrix_index=count++;
        matrix_index_to_cell_index(matrix_index)=cell_index;}
    for(int i=0;i<row_counts.m;i++){
        int boundary=0;
        for(int j=0;j<TV::m;j++) if(matrix_index_to_cell_index(i)(j)==0 || matrix_index_to_cell_index(i)(j)==mac_grid.counts(j)-1) boundary++;
        row_counts(i)=(2*TV::m+1)-boundary;}
    A.Set_Row_Lengths(row_counts);
    TV one_over_dx2=Inverse(mac_grid.dX*mac_grid.dX);
    T default_row_sum=-2*one_over_dx2.Sum_Abs();
    TV_INT grid_counts=mac_grid.counts;
    for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        T row_sum=default_row_sum;
        int matrix_index=cell_index_to_matrix_index(cell_index);
        for(int axis=0;axis<TV::m;axis++){TV_INT offset;offset[axis]=1;
            if(elliptic_solver->psi_N.Component(axis)(cell_index)) row_sum+=one_over_dx2[axis];
            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),one_over_dx2[axis]);
            if(elliptic_solver->psi_N.Component(axis)(cell_index+offset)) row_sum+=one_over_dx2[axis];
            else A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),one_over_dx2[axis]);}
        A.Set_Element(matrix_index,matrix_index,row_sum);}
}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class TV> void FAST_PROJECTION_DYNAMICS_UNIFORM<TV>::
Make_Divergence_Free_Fast(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    Compute_Divergence(FACE_LOOKUP_UNIFORM<TV>(face_velocities),elliptic_solver);
    for(CELL_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        int matrix_index=cell_index_to_matrix_index(cell_index);
        b(matrix_index)=elliptic_solver->f(cell_index);}
    elliptic_solver->Solve_Subregion(matrix_index_to_cell_index,A,b);
    Apply_Pressure(face_velocities,dt,time);
    A.Negate();
}
//#####################################################################
namespace PhysBAM{
template class FAST_PROJECTION_DYNAMICS_UNIFORM<VECTOR<float,1> >;
template class FAST_PROJECTION_DYNAMICS_UNIFORM<VECTOR<float,2> >;
template class FAST_PROJECTION_DYNAMICS_UNIFORM<VECTOR<float,3> >;
template class FAST_PROJECTION_DYNAMICS_UNIFORM<VECTOR<double,1> >;
template class FAST_PROJECTION_DYNAMICS_UNIFORM<VECTOR<double,2> >;
template class FAST_PROJECTION_DYNAMICS_UNIFORM<VECTOR<double,3> >;
}
