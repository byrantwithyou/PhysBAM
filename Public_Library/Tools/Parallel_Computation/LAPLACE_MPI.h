//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_MPI
//#####################################################################
#ifndef __LAPLACE_MPI__
#define __LAPLACE_MPI__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Parallel_Computation/MPI_GRID.h>
#include <Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <Tools/Parallel_Computation/PCG_SPARSE_THREADED.h>
#include <Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T> class PCG_SPARSE;
template<class T> class MPI_UNIFORM_GRID;
template<class TV> struct GRID_ARRAYS_POLICY;
class SPARSE_MATRIX_PARTITION;
template<class TV> class LAPLACE_UNIFORM;

template<class TV>
class LAPLACE_MPI:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef TV_INT T_INDEX;typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename ARRAY<T,TV_INT>::template REBIND<int>::TYPE T_ARRAYS_INT;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef GRID<TV> T_PARALLEL_GRID;
public:
    MPI_UNIFORM_GRID<TV>*& mpi_grid;
    const GRID<TV>& local_grid;
    PCG_SPARSE<T>& local_pcg;
    PCG_SPARSE_THREADED<TV>* local_pcg_threaded;
    int& number_of_regions;
    int number_of_global_regions;
    T_ARRAYS_INT& filled_region_colors;
    ARRAY<bool>& filled_region_touches_dirichlet;
    bool& solve_neumann_regions;
    T_FACE_ARRAYS_BOOL& psi_N;
    ARRAY<ARRAY<int> > filled_region_ranks;
    ARRAY<SPARSE_MATRIX_PARTITION> partitions;
    ARRAY<MPI::Group>* groups;
    ARRAY<MPI::Intracomm>* communicators;

    LAPLACE_MPI(LAPLACE_UNIFORM<TV>& laplace);
    virtual ~LAPLACE_MPI();


//#####################################################################
    void Synchronize_Solution_Regions();
    void Update_Solution_Regions_For_Solid_Fluid_Coupling(const MPI_UNIFORM_GRID<TV>& mpi_grid);
    virtual void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<T_INDEX> >& matrix_index_to_cell_index_array)=0;
    int Get_Total_Number_Of_Threads(const int input,const int color);
    void Solve_Threaded(RANGE<TV_INT>& domain,const ARRAY<int,TV_INT>& domain_index,ARRAY<INTERVAL<int> >& interior_indices,ARRAY<ARRAY<INTERVAL<int> > >& ghost_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,ARRAY<KRYLOV_VECTOR_BASE<T>*>& vectors,const T tolerance,const int color,const int multi_proc_mode);
    void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,ARRAY<KRYLOV_VECTOR_BASE<T>*>& vectors,const T tolerance,const int color);
    void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,const T tolerance,const int color,const ARRAY<VECTOR<int,2> >& global_column_index_boundaries);
    static bool& Use_Parallel_Solve();
//#####################################################################
};
}
#endif
