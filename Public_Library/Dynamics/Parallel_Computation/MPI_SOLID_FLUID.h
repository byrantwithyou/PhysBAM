//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_SOLID_FLUID
//#####################################################################
#ifndef __MPI_SOLID_FLUID__
#define __MPI_SOLID_FLUID__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

class MPI_PACKAGE;
template<class TV> class FLUID_SYSTEM_MPI;
template<class TV> class SOLID_SYSTEM_MPI;
template<class TV> class GENERALIZED_VELOCITY;
template<class TV> class SOLID_BODY_COLLECTION;
class SPARSE_MATRIX_PARTITION;

template<class TV>
class MPI_SOLID_FLUID:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    int rank; // global rank
    int number_of_processes;
    int number_of_solid_processes;
    int solid_node;
    MPI::Intracomm* comm;
    MPI::Group *group,*solid_group,*fluid_group;
    ARRAY<int> solid_ranks,fluid_ranks;
private:
    mutable int current_tag;
public:

    MPI_SOLID_FLUID();
    ~MPI_SOLID_FLUID();

    int Number_Of_Processors() const
    {return number_of_processes;}

    int Get_Unique_Tag() const
    {return current_tag=max(32,(current_tag+1)&((1<<15)-1));}

public:
    void Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<TV>& solid_body_collection) const;
    bool Fluid_Node() const;
    bool Solid_Node() const;
    void Create_Fluid_Comm_For_Solid_Nodes() const;
    template<class T2> void Reduce_Add(const T2& input,T2& output) const;
    T Reduce_Min(const T local_value) const;
    T Reduce_Max(const T local_value) const;
    void Parallel_Solve_Fluid_Part(FLUID_SYSTEM_MPI<TV>& fluid_system,KRYLOV_VECTOR_WRAPPER<T,ARRAY<ARRAY<T> > >& x_array,
        KRYLOV_VECTOR_WRAPPER<T,ARRAY<ARRAY<T> > >& b_array,ARRAY<KRYLOV_VECTOR_BASE<T>*>& vectors,
        const int min_iterations,const int max_iterations,const T tolerance,const bool recompute_preconditioner,
        ARRAY<MPI::Intracomm>* fluid_comm,ARRAY<SPARSE_MATRIX_PARTITION>* partitions);
    void Parallel_Solve_Solid_Part(SOLID_SYSTEM_MPI<TV>& solid_system,GENERALIZED_VELOCITY<TV>& x_array,GENERALIZED_VELOCITY<TV>& b_array,
        ARRAY<KRYLOV_VECTOR_BASE<T>*>& vectors,const int min_iterations,const int max_iterations,const T tolerance);
    void Distribute_Lists_From_Solid_Node(GENERALIZED_VELOCITY<TV>& F) const;
    void Exchange_Coupled_Deformable_Particle_List(ARRAY<int>* fluid_list,ARRAY<ARRAY<int> >* results);
    void Aggregate_Lists_To_Solid_Node(GENERALIZED_VELOCITY<TV>& F);
//#####################################################################
};
}
#endif
