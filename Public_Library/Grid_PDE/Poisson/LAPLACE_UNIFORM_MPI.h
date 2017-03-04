//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_UNIFORM_MPI
//#####################################################################
#ifndef __LAPLACE_UNIFORM_MPI__
#define __LAPLACE_UNIFORM_MPI__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/RANGE.h>
#include <Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Grid_PDE/Poisson/LAPLACE_MPI.h>
namespace PhysBAM{

class GRAPH;
template<class TV> class GRID;
template<class TV> class LAPLACE_UNIFORM;

template<class TV>
class LAPLACE_UNIFORM_MPI:public LAPLACE_MPI<TV>
{
    typedef typename TV::SCALAR T;
    typedef MPI_UNIFORM_GRID<TV> T_MPI_GRID;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
public:
    typedef LAPLACE_MPI<TV> BASE;
    using BASE::mpi_grid;using BASE::local_grid;using BASE::filled_region_ranks;using BASE::partitions;using BASE::number_of_regions;using BASE::solve_neumann_regions;using BASE::psi_N;
    using BASE::filled_region_touches_dirichlet;

    bool sum_int_needs_init;
    int sum_int;

    LAPLACE_UNIFORM_MPI(LAPLACE_UNIFORM<TV>& laplace)
        :LAPLACE_MPI<TV>(laplace),sum_int_needs_init(false),sum_int(0)
    {
    }

    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array) override
    {Find_Matrix_Indices(filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array,TV());}

    int Global_Sum(int input)
    {
        return input;
    }

//#####################################################################
private:
    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,1>&);
    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,2>&);
    void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,const VECTOR<T,3>&);
    void Find_Matrix_Indices_In_Region(const int region_index,const RANGE<TV_INT>& region,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,
        ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array);
    void Find_Boundary_Indices_In_Region(const int side,const RANGE<TV_INT>& region,T_ARRAYS_INT& cell_index_to_matrix_index);
//#####################################################################
};
}
#endif
