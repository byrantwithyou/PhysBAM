//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_GRID
//#####################################################################
#ifndef __MPI_GRID__
#define __MPI_GRID__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Vectors/VECTOR_2D.h>
namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

class MPI_PACKAGE;
template<class TV> class RANGE;
template<class TV> class GRID;

template<class TV>
class MPI_GRID:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
public:
    GRID<TV>& local_grid;
    int number_of_ghost_cells;
    GRID<TV> global_grid;
    int rank;
    int number_of_processes;
    MPI::Intracomm* comm;
    MPI::Group* group;
    GRID<TV> process_grid;
    T_ARRAYS_INT process_ranks;
    TV_INT coordinates;
    ARRAY<TV_INT> all_coordinates;
    ARRAY<int> all_neighbor_ranks,side_neighbor_ranks; // all_neighbor_ranks includes all 1 ring neighbors
    ARRAY<TV_INT> all_neighbor_directions,side_neighbor_directions;
    TV_INT local_to_global_offset;
    T_ARRAYS_INT local_cell_index_to_global_column_index_map;
    ARRAY<ARRAY<int> > boundaries;
private:
    mutable int current_tag;
public:
    VECTOR<bool,TV::m> periodic;
    bool ignore_boundary_faces;

    MPI_GRID(GRID<TV>& local_grid_input,const int number_of_ghost_cells_input,const bool skip_initialization=false,const TV_INT& processes_per_dimension=TV_INT(),
        const VECTOR<bool,TV::m>& periodic_input=(VECTOR<bool,TV::m>()),MPI::Group* group_input=0);
    ~MPI_GRID();

    int Number_Of_Processors() const
    {return number_of_processes;}

    int Get_Unique_Tag() const
    {return current_tag=max(32,(current_tag+1)&((1<<15)-1));}

    int Get_Send_Tag(const TV_INT& direction) const
    {STATIC_ASSERT(TV_INT::m<=3);int tag=0;
    for(int i=0;i<direction.m;i++){assert(abs(direction[i])<=1);tag=3*tag+direction[i]+1;}
    return tag;}

    int Get_Recv_Tag(const TV_INT& direction) const
    {return Get_Send_Tag(-direction);}

protected:
    TV Wrap_Offset(const TV_INT& direction) const // offset to add to translate into space of adjacent processor
    {TV offset;TV_INT neighbor_coordinates=coordinates+direction;
    for(int axis=0;axis<offset.m;axis++)if(periodic[axis] && (neighbor_coordinates[axis]<1 || neighbor_coordinates[axis]>process_grid.counts[axis]))
        offset[axis]=-direction[axis]*global_grid.domain.Edge_Lengths()[axis];
    return offset;}
public:

//#####################################################################
    void Initialize(VECTOR<VECTOR<bool,2>,TV::m>& domain_walls);
    bool Neighbor(const int axis,const int axis_side) const;
    void Split_Grid(const TV_INT& processes_per_dimension);
    GRID<TV> Restrict_Grid(const TV_INT& coordinates) const;
    void Synchronize_Dt(T& dt) const;
    void Synchronize_J_Bounds(int& jmin,int& jmax) const;
    void Sync_Common_Face_Weights_To(ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >,FACE_INDEX<TV::dimension> >& weights_to,ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >,FACE_INDEX<TV::dimension> >& weights_from,const int ghost_cells);
    void Sync_Common_Face_Weights_From(ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >,FACE_INDEX<TV::dimension> >& weights_to,ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >,FACE_INDEX<TV::dimension> >& weights_from,const int ghost_cells);
    void Sync_Common_Cell_Weights_To(ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& weights_to,ARRAY<ARRAY<PAIR<TV_INT,int> >,TV_INT>& weights_from,const int ghost_cells);
    void Sync_Common_Cell_Weights_From(ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& weights_to,ARRAY<ARRAY<PAIR<TV_INT,int> >,TV_INT>& weights_from,const int ghost_cells);
    template<class T_MPI_GRID,class T2> void Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,ARRAYS_ND_BASE<T2,VECTOR<int,TV::dimension> >& data,const int bandwidth,const bool include_corners=true) const;
    template<class T_MPI_GRID,class T2> void Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,const GRID<TV>& local_grid,ARRAYS_ND_BASE<T2,VECTOR<int,TV::dimension> >& data,const int bandwidth,
        const bool include_corners=true) const;
    template<class T_MPI_GRID,class T2,class T_ARRAYS,class INDEX> void Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,ARRAY_BASE<T2,T_ARRAYS,INDEX>& data,const int bandwidth,const bool include_corners=true) const;
    template<class T_MPI_GRID,class T2,class T_ARRAYS,class INDEX> void Exchange_Boundary_Cell_Data(const T_MPI_GRID& mpi_grid,const GRID<TV>& local_grid,ARRAY_BASE<T2,T_ARRAYS,INDEX>& data,const int bandwidth,
        const bool include_corners=true) const;
    template<class T_MPI_GRID,class T_FACE_ARRAYS> void Exchange_Boundary_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data,const int bandwidth) const;
    template<class T_MPI_GRID,class T_FACE_ARRAYS> void Average_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data) const;
    template<class T_MPI_GRID,class T_FACE_ARRAYS> void Copy_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data) const;
    template<class T_MPI_GRID,class T_FACE_ARRAYS> void Assert_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_ARRAYS& data,const T tolerance=0) const;
    template<class T_MPI_GRID,class T_FACE_SCALAR> void Union_Common_Face_Data(const T_MPI_GRID& mpi_grid,T_FACE_SCALAR& data) const;
    void Find_Boundary_Regions(ARRAY<RANGE<TV_INT> >& regions,const RANGE<TV_INT>& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
        const bool include_ghost_regions=true) const;
    void Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,1> > >& regions,const RANGE<VECTOR<int,1> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
        const bool include_ghost_regions,const GRID<VECTOR<T,1> >& local_grid) const;
    void Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,2> > >& regions,const RANGE<VECTOR<int,2> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
        const bool include_ghost_regions,const GRID<VECTOR<T,2> >& local_grid) const;
    void Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,3> > >& regions,const RANGE<VECTOR<int,3> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
        const bool include_ghost_regions,const GRID<VECTOR<T,3> >& local_grid) const;
    RANGE<TV_INT> Find_Region_Box(const int processor,const RANGE<TV_INT>& sentinels,const int band) const;
    template<class T2> void Reduce_Add(const T2& input,T2& output) const;
    template<class T2> T2 Reduce_Add(const T2& local_value) const;

protected:
    VECTOR<int,1> Split_Grid(const GRID<VECTOR<T,1> >& global_grid,const VECTOR<int,1>& processes_per_dimension);
    VECTOR<int,2> Split_Grid(const GRID<VECTOR<T,2> >& global_grid,const VECTOR<int,2>& processes_per_dimension);
    VECTOR<int,3> Split_Grid(const GRID<VECTOR<T,3> >& global_grid,const VECTOR<int,3>& processes_per_dimension);
    static void Split_Dimension(const int m,const int processes,ARRAY<int>& boundaries);
    void Initialize_Communicator(const bool manual,MPI::Group* group);
//#####################################################################
};
}
#endif
