//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PCG_SPARSE_MPI
//#####################################################################
#ifndef __PCG_SPARSE_MPI__
#define __PCG_SPARSE_MPI__

#ifdef USE_MPI

#include <Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Parallel_Computation/THREADED_UNIFORM_GRID.h>
namespace PhysBAM{

class SPARSE_MATRIX_PARTITION;

template<class TV>
class PCG_SPARSE_MPI:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    PCG_SPARSE<T>& pcg;
    MPI::Intracomm& comm;
    THREADED_UNIFORM_GRID<TV>* thread_grid;
    SPARSE_MATRIX_PARTITION& partition;
    ARRAY<MPI::Datatype> boundary_datatypes,ghost_datatypes;
    ARRAY<ARRAY<int> > columns_to_send;
    ARRAY<ARRAY<int> > columns_to_receive;

    PCG_SPARSE_MPI(PCG_SPARSE<T>& pcg_input,MPI::Intracomm& comm_input,SPARSE_MATRIX_PARTITION& partition_input)
        :pcg(pcg_input),comm(comm_input),partition(partition_input)
    {}

    ~PCG_SPARSE_MPI()
    {MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes);MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes);}

    template<class TYPE> TYPE Global_Sum(const TYPE& input)
    {TYPE output;MPI_UTILITIES::Reduce(input,output,MPI::SUM,comm);return output;}

    template<class TYPE> TYPE Global_Max(const TYPE& input)
    {TYPE output;MPI_UTILITIES::Reduce(input,output,MPI::MAX,comm);return output;}

//#####################################################################
    void Fill_Ghost_Cells(ARRAY<T>& v);
    void Serial_Solve(SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,ARRAY<T>& q,ARRAY<T>& s,ARRAY<T>& r,ARRAY<T>& k,ARRAY<T>& z,const int tag,const T tolerance=1e-7);
    void Parallel_Solve(SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,const T tolerance=1e-7,const bool recompute_preconditioner=true);
    void Parallel_Solve(SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& x_local,ARRAY<T>& b_local,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries,
        const T tolerance=1e-7,const bool recompute_preconditioner=true);
    void Find_Ghost_Regions(SPARSE_MATRIX_FLAT_MXN<T>& A,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries);
    void Find_Ghost_Regions_Threaded(SPARSE_MATRIX_FLAT_MXN<T>& A,const ARRAY<VECTOR<int,2> >& proc_column_index_boundaries);
    void Fill_Ghost_Cells_Far(ARRAY<T>& x);
    void Fill_Ghost_Cells_Threaded(ARRAY<T>& x);    
    void Initialize_Datatypes();
//#####################################################################
};
}
#endif
#endif
