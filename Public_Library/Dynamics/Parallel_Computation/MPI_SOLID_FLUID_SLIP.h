//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_SOLID_FLUID_SLIP
//#####################################################################
#ifndef __MPI_SOLID_FLUID_SLIP__
#define __MPI_SOLID_FLUID_SLIP__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS_BINARY_UNIFORM.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Parallel_Computation/MPI_GRID_POLICY.h>
#include <Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

class MPI_PACKAGE;
template<class TV> class FLUID_SYSTEM_MPI_SLIP;
template<class TV> class SOLID_SYSTEM_MPI_SLIP;
template<class TV> class GENERALIZED_VELOCITY;
template<class TV> class SOLID_BODY_COLLECTION;
template<class T_GRID> class POISSON_COLLIDABLE_UNIFORM;
template<class TV> class GRID;

template<class TV>
class MPI_SOLID_FLUID_SLIP:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename GRID<TV>::VECTOR_INT TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef ARRAY<T,SIDED_FACE_INDEX<TV::dimension> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename MPI_GRID_POLICY<GRID<TV> >::MPI_GRID T_MPI_GRID;
    typedef typename MPI_GRID_POLICY<GRID<TV> >::PARALLEL_GRID T_PARALLEL_GRID;
public:
    int rank; // global rank
    int number_of_processes;
    int number_of_solid_processes;
    int solid_node;
    MPI::Intracomm* comm;
    MPI::Intracomm* fluid_comm;
    MPI::Group *group,*solid_group,*fluid_group;
    ARRAY<int> solid_ranks,fluid_ranks;
    SPARSE_MATRIX_PARTITION partition;
    T_MPI_GRID* mpi_grid;
    POISSON_COLLIDABLE_UNIFORM<GRID<TV> >* poisson;
private:
    mutable int current_tag;
public:

    MPI_SOLID_FLUID_SLIP();
    ~MPI_SOLID_FLUID_SLIP();

    int Number_Of_Processors() const
    {return number_of_processes;}

    int Get_Unique_Tag() const
    {current_tag=max(32,(current_tag+1)&((1<<15)-1));return current_tag;}

    void Set_Poisson(POISSON_COLLIDABLE_UNIFORM<GRID<TV> >& poisson_input)
    {poisson=&poisson_input;}

    void Set_MPI_Grid(T_MPI_GRID& mpi_grid_input);

public:
    void Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<TV>& solid_body_collection) const;
    bool Fluid_Node() const;
    bool Solid_Node() const;
    void Create_Fluid_Comm_For_Solid_Nodes() const;
    template<class T2> void Reduce_Add(const T2& input,T2& output) const;
    T Reduce_Min(const T local_value) const;
    void Parallel_Solve_Fluid_Part(FLUID_SYSTEM_MPI_SLIP<TV>& fluid_system,ARRAY<T>& x_array,ARRAY<T>& b_array,ARRAY<T>& p_array,ARRAY<T>& ap_array,
        ARRAY<T>& ar_array,ARRAY<T>& r_array,ARRAY<T>& z_array,ARRAY<T>& zaq_array,const int min_iterations,const int max_iterations,
        const T tolerance,const bool recompute_preconditioner,const bool leakproof_solve);
    void Parallel_Solve_Solid_Part(SOLID_SYSTEM_MPI_SLIP<TV>& solid_system,GENERALIZED_VELOCITY<TV>& x_array,GENERALIZED_VELOCITY<TV>& b_array,GENERALIZED_VELOCITY<TV>& p_array,GENERALIZED_VELOCITY<TV>& ap_array,
        GENERALIZED_VELOCITY<TV>& ar_array,GENERALIZED_VELOCITY<TV>& r_array,GENERALIZED_VELOCITY<TV>& z_array,GENERALIZED_VELOCITY<TV>& zaq_array,const int min_iterations,const int max_iterations,const T tolerance);
    void Exchange_Coupled_Deformable_Particle_List(ARRAY<int>* fluid_list,ARRAY<ARRAY<int> >* results);
    void Find_Matrix_Indices(const GRID<TV>& local_grid,const ARRAY<bool,TV_INT>& valid_divergence_cells,T_ARRAYS_INT& cell_index_to_matrix_index,T_FACE_ARRAYS_INT& face_ghost_cell_index,ARRAY<int>& face_ghost_cell_index_map,T_FACE_ARRAYS_INT& face_lambdas,INTERVAL<int>& divergence_indices,int &cell_count);
    //void Find_Matrix_Indices_In_Region(const GRID<TV>& local_grid,const ARRAY<bool,TV_INT>& valid_divergence_cells,const int region_index,const RANGE<TV_INT>& region,T_ARRAYS_INT& cell_index_to_matrix_index,
    //    T_FACE_ARRAYS_INT& face_ghost_cell_index,ARRAY<int>& face_ghost_cell_index_map,T_FACE_ARRAYS_INT& face_lambdas,INTERVAL<int>& divergence_indices,int& cell_count);
    void Find_Boundary_Indices_In_Region(const GRID<TV>& local_grid,const ARRAY<bool,TV_INT>& valid_divergence_cells,const int domain_side,const RANGE<TV_INT>& region,const T_ARRAYS_INT& cell_index_to_matrix_index,const T_FACE_ARRAYS_INT& face_ghost_cell_index,const T_FACE_ARRAYS_INT& face_lambdas);
//#####################################################################
};
}
#endif
