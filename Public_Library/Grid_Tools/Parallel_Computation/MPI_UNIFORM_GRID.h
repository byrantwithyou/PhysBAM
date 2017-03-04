//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_UNIFORM_GRID
//#####################################################################
#ifndef __MPI_UNIFORM_GRID__
#define __MPI_UNIFORM_GRID__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Grid_Tools/Parallel_Computation/MPI_GRID.h>
namespace PhysBAM{

template<class TV>
class MPI_UNIFORM_GRID:public MPI_GRID<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef ARRAY<RANGE<TV_INT> ,TV_INT> T_ARRAYS_BOX_INT;
public:
    typedef GRID<TV> GRID_T;

    typedef MPI_GRID<TV> BASE;
    using BASE::Exchange_Boundary_Cell_Data;using BASE::local_grid;using BASE::global_grid;using BASE::comm;using BASE::all_neighbor_ranks;using BASE::side_neighbor_ranks;
    using BASE::Get_Unique_Tag;using BASE::Get_Send_Tag;using BASE::Get_Recv_Tag;using BASE::Wrap_Offset;using BASE::all_neighbor_directions;using BASE::all_coordinates;
    using BASE::Restrict_Grid;

    MPI_UNIFORM_GRID(GRID<TV>& local_grid_input,const int number_of_ghost_cells_input,const bool skip_initialization=false,const TV_INT& processes_per_dimension=TV_INT(),
        const VECTOR<bool,TV::m>& periodic_input=(VECTOR<bool,TV::m>()),MPI::Group* group_input=0);
    
    template<class T_ARRAYS> void Exchange_Boundary_Cell_Data(T_ARRAYS& data,const int bandwidth,const bool include_corners=true) const
    {MPI_GRID<TV>::Exchange_Boundary_Cell_Data(*this,data,bandwidth,include_corners);}
    
    template<class T_ARRAYS> void Union_Common_Face_Data(T_ARRAYS& data) const
    {MPI_GRID<TV>::Union_Common_Face_Data(*this,data);}

    template<class T_FACE_ARRAYS_2>
    void Exchange_Boundary_Face_Data(T_FACE_ARRAYS_2& data,const int bandwidth) const
    {MPI_GRID<TV>::Exchange_Boundary_Face_Data(*this,data,bandwidth);}

    template<class T_FACE_ARRAYS_2> void Average_Common_Face_Data(T_FACE_ARRAYS_2& data) const
    {MPI_GRID<TV>::Average_Common_Face_Data(*this,data);}

    template<class T_FACE_ARRAYS_2> void Copy_Common_Face_Data(T_FACE_ARRAYS_2& data) const
    {MPI_GRID<TV>::Copy_Common_Face_Data(*this,data);}

    template<class T_FACE_ARRAYS_2> void Assert_Common_Face_Data(T_FACE_ARRAYS_2& data,const T tolerance=0) const
    {MPI_GRID<TV>::Assert_Common_Face_Data(*this,data,tolerance);}

    RANGE<TV_INT> Face_Sentinels(const int axis) const
    {return RANGE<TV_INT>(TV_INT(),TV_INT::Axis_Vector(axis));}

    RANGE<TV_INT> Parallel_Face_Sentinels(const int axis) const
    {return Face_Sentinels(axis);}
    
    void Synchronize_Dt(T& dt) const
    {MPI_GRID<TV>::Synchronize_Dt(dt);}
    
    void Initialize(VECTOR<VECTOR<bool,2>,TV::m>& domain_walls)
    {BASE::Initialize(domain_walls);}
    
    bool Neighbor(const int axis,const int axis_side) const
    {return BASE::Neighbor(axis,axis_side);}

//#####################################################################
    GRID<TV> Get_Non_Overlapping_Face_Grid(const int axis) const;
    template<class T_ARRAYS> bool Gather_Cell_Data(const T_ARRAYS& local_data,T_ARRAYS& global_data) const;
    template<class T_ARRAYS> void Scatter_Cell_Data(const T_ARRAYS& global_data,T_ARRAYS& local_data) const;
    template<class T2> MPI_PACKAGE Package_Cell_Data(ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& data,const RANGE<TV_INT>& region) const;
    template<class T_FACE_ARRAYS1> MPI_PACKAGE Package_Face_Data(T_FACE_ARRAYS1& data,const ARRAY<RANGE<TV_INT> >& region) const;
    template<class T_FACE_ARRAYS1> MPI_PACKAGE Package_Common_Face_Data(T_FACE_ARRAYS1& data,const int axis,const RANGE<TV_INT>& region) const;
//#####################################################################
};
}
#endif
