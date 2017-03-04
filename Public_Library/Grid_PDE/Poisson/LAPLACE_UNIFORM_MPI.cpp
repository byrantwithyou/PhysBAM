//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Frank Losasso, Michael Lentine, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_PDE/Poisson/LAPLACE_UNIFORM.h>
#include <Grid_PDE/Poisson/LAPLACE_UNIFORM_MPI.h>
#ifdef USE_MPI
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Data_Structures/GRAPH.h>
#include <Core/Data_Structures/UNION_FIND.h>
#include <Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Grid_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
#endif
using namespace PhysBAM;

#ifdef USE_MPI

//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class TV> void LAPLACE_UNIFORM_MPI<TV>::
Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,1>&)
{
    assert(local_grid.Is_MAC_Grid());
    int m=local_grid.counts.x;
    Find_Matrix_Indices_In_Region(0,RANGE<VECTOR<int,1> >(1,m),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(1,RANGE<VECTOR<int,1> >(0,0),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(2,RANGE<VECTOR<int,1> >(m+1,m+1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Boundary_Indices_In_Region(1,RANGE<VECTOR<int,1> >(1,1),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(2,RANGE<VECTOR<int,1> >(m,m),cell_index_to_matrix_index);
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class TV> void LAPLACE_UNIFORM_MPI<TV>::
Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,2>&)
{
    assert(local_grid.Is_MAC_Grid());
    int m=local_grid.counts.x,n=local_grid.counts.y;
    Find_Matrix_Indices_In_Region(0,RANGE<VECTOR<int,2> >(1,m,1,n),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(1,RANGE<VECTOR<int,2> >(0,0,1,n),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(2,RANGE<VECTOR<int,2> >(m+1,m+1,1,n),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(3,RANGE<VECTOR<int,2> >(1,m,0,0),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(4,RANGE<VECTOR<int,2> >(1,m,n+1,n+1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Boundary_Indices_In_Region(1,RANGE<VECTOR<int,2> >(1,1,1,n),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(2,RANGE<VECTOR<int,2> >(m,m,1,n),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(3,RANGE<VECTOR<int,2> >(1,m,1,1),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(4,RANGE<VECTOR<int,2> >(1,m,n,n),cell_index_to_matrix_index);
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class TV> void LAPLACE_UNIFORM_MPI<TV>::
Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,3>&)
{
    assert(local_grid.Is_MAC_Grid());
    int m=local_grid.counts.x,n=local_grid.counts.y,mn=local_grid.counts.z;
    Find_Matrix_Indices_In_Region(0,RANGE<VECTOR<int,3> >(1,m,1,n,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(1,RANGE<VECTOR<int,3> >(0,0,1,n,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(2,RANGE<VECTOR<int,3> >(m+1,m+1,1,n,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(3,RANGE<VECTOR<int,3> >(1,m,0,0,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(4,RANGE<VECTOR<int,3> >(1,m,n+1,n+1,1,mn),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(5,RANGE<VECTOR<int,3> >(1,m,1,n,0,0),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Matrix_Indices_In_Region(6,RANGE<VECTOR<int,3> >(1,m,1,n,mn+1,mn+1),filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    Find_Boundary_Indices_In_Region(1,RANGE<VECTOR<int,3> >(1,1,1,n,1,mn),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(2,RANGE<VECTOR<int,3> >(m,m,1,n,1,mn),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(3,RANGE<VECTOR<int,3> >(1,m,1,1,1,mn),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(4,RANGE<VECTOR<int,3> >(1,m,n,n,1,mn),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(5,RANGE<VECTOR<int,3> >(1,m,1,n,1,1),cell_index_to_matrix_index);
    Find_Boundary_Indices_In_Region(6,RANGE<VECTOR<int,3> >(1,m,1,n,mn,mn),cell_index_to_matrix_index);
}
//#####################################################################
// Function Find_Matrix_Indices_In_Region
//#####################################################################
template<class TV> void LAPLACE_UNIFORM_MPI<TV>::
Find_Matrix_Indices_In_Region(const int region_index,const RANGE<TV_INT>& region,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,
    ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array)
{
    if(region_index>=0) for(int color=0;color<filled_region_ranks.m;color++)partitions(color).ghost_indices(region_index).min_corner=filled_region_cell_count(color);
    else for(int color=0;color<filled_region_ranks.m;color++)partitions(color).interior_indices.min_corner=filled_region_cell_count(color);
    for(CELL_ITERATOR<TV> iterator(local_grid,region);iterator.Valid();iterator.Next()){TV_INT c=iterator.Cell_Index();
        int color=filled_region_colors(c);if(color<1 || (!filled_region_touches_dirichlet(color)&&!solve_neumann_regions)) continue;
        int new_index=filled_region_cell_count(color)++;cell_index_to_matrix_index(c)=new_index;
        matrix_index_to_cell_index_array(color)(new_index)=c;}
    if(region_index>=0) for(int color=0;color<filled_region_ranks.m;color++)partitions(color).ghost_indices(region_index).max_corner=filled_region_cell_count(color);
    else for(int color=0;color<filled_region_ranks.m;color++)partitions(color).interior_indices.max_corner=filled_region_cell_count(color);
}
//#####################################################################
// Function Find_Boundary_Indices_In_Region
//#####################################################################
template<class TV> void LAPLACE_UNIFORM_MPI<TV>::
Find_Boundary_Indices_In_Region(const int side,const RANGE<TV_INT>& region,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    int axis=side/2,cell_side=side&1;
    RANGE<TV_INT> face_region=region;if(cell_side) face_region+=TV_INT::Axis_Vector(axis);
    // count boundary indices
    ARRAY<int> counts(partitions.m);
    TV_INT face_offset=cell_side?TV_INT::Axis_Vector(axis):TV_INT();
    for(CELL_ITERATOR<TV> iterator(local_grid,region);iterator.Valid();iterator.Next())if(!psi_N.Component(axis)(iterator.Cell_Index()+face_offset)){
        int color=filled_region_colors(iterator.Cell_Index());
        if(counts.Valid_Index(color)) counts(color)++;}
    // fill boundary indices
    for(int color=0;color<partitions.m;color++)partitions(color).boundary_indices(side).Resize(counts(color));
    counts.Fill(0);
    for(CELL_ITERATOR<TV> iterator(local_grid,region);iterator.Valid();iterator.Next())if(!psi_N.Component(axis)(iterator.Cell_Index()+face_offset)){
        int color=filled_region_colors(iterator.Cell_Index());
        if(counts.Valid_Index(color)) partitions(color).boundary_indices(side)(counts(color)++)=cell_index_to_matrix_index(iterator.Cell_Index());}
}
//#####################################################################

#else

//#####################################################################
template<class TV> void LAPLACE_UNIFORM_MPI<TV>::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,const VECTOR<T,1>&){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void LAPLACE_UNIFORM_MPI<TV>::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,const VECTOR<T,2>&){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void LAPLACE_UNIFORM_MPI<TV>::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,const VECTOR<T,3>&){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################

#endif

//#####################################################################
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER(TV) \
    template void LAPLACE_UNIFORM_MPI<TV>::Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >&,T_ARRAYS_INT&,ARRAY<ARRAY<TV_INT> >&,const TV&);
template LAPLACE_UNIFORM_MPI<VECTOR<float,1> >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(VECTOR<float,1>) >&);
template LAPLACE_UNIFORM_MPI<VECTOR<float,2> >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(VECTOR<float,2>) >&);
template LAPLACE_UNIFORM_MPI<VECTOR<float,3> >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(VECTOR<float,3>) >&);
INSTANTIATION_HELPER(P(VECTOR<float,1>));
INSTANTIATION_HELPER(P(VECTOR<float,2>));
INSTANTIATION_HELPER(P(VECTOR<float,3>));
template LAPLACE_UNIFORM_MPI<VECTOR<double,1> >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(VECTOR<double,1>) >&);
template LAPLACE_UNIFORM_MPI<VECTOR<double,2> >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(VECTOR<double,2>) >&);
template LAPLACE_UNIFORM_MPI<VECTOR<double,3> >::LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<P(VECTOR<double,3>) >&);
INSTANTIATION_HELPER(P(VECTOR<double,1>));
INSTANTIATION_HELPER(P(VECTOR<double,2>));
INSTANTIATION_HELPER(P(VECTOR<double,3>));
