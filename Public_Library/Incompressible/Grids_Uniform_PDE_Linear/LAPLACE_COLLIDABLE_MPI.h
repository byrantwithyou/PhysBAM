//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_COLLIDABLE_MPI
//#####################################################################
#ifndef __LAPLACE_COLLIDABLE_MPI__
#define __LAPLACE_COLLIDABLE_MPI__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Parallel_Computation/MPI_GRID.h>
#include <Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/LAPLACE_COLLIDABLE_POLICY.h>
namespace PhysBAM{

template<class T> class PCG_SPARSE;
template<class TV> struct GRID_ARRAYS_POLICY;
template<class T> class MPI_UNIFORM_GRID;

template<class TV>
class LAPLACE_COLLIDABLE_MPI:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;typedef MPI_UNIFORM_GRID<TV> T_MPI_GRID;
    typedef VECTOR<int,TV::m> TV_INT;typedef typename LAPLACE_COLLIDABLE_POLICY<TV>::LAPLACE T_LAPLACE_COLLIDABLE;
    typedef TV_INT T_INDEX;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef GRID<TV> T_PARALLEL_GRID;
public:
    T_MPI_GRID*& mpi_grid;
    const GRID<TV>& local_grid;
    PCG_SPARSE<T>& local_pcg;
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

    LAPLACE_COLLIDABLE_MPI(T_LAPLACE_COLLIDABLE& laplace);
    virtual ~LAPLACE_COLLIDABLE_MPI();

//#####################################################################
    void Synchronize_Solution_Regions();
    void Update_Solution_Regions_For_Solid_Fluid_Coupling(const T_MPI_GRID& mpi_grid);
    virtual void Find_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index,ARRAY<ARRAY<T_INDEX> >& matrix_index_to_cell_index_array)=0;
    void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,ARRAY<T>& q,ARRAY<T>& s,ARRAY<T>& r,ARRAY<T>& k,ARRAY<T>& z,const T tolerance,const int color);
    void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,ARRAY<T>& x,ARRAY<T>& b,const T tolerance,const int color,const ARRAY<VECTOR<int,2> >& global_column_index_boundaries);
    static bool& Use_Parallel_Solve();
//#####################################################################
};
}
#endif
