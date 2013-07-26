//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LOCAL_GRID
//#####################################################################
#ifndef __LOCAL_GRID__
#define __LOCAL_GRID__

#include <Tools/Log/LOG.h>
#include <Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
namespace PhysBAM{

template<class TV>
class LOCAL_GRID
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
public:
    typedef int HAS_UNTYPED_READ_WRITE;

    GRID<TV> grid;
    MPI_UNIFORM_GRID<TV> mpi_grid;
    const GRID<TV>& global_grid;
    TV_INT offset;
    ARRAY<bool> neighbor_overlaps;

    LOCAL_GRID(const GRID<TV>& global_grid_input);
    LOCAL_GRID(const GRID<TV>& global_grid_input,const GRID<TV>& local_grid_input);
    ~LOCAL_GRID();

    RANGE<TV_INT> Interior_Region(const RANGE<TV_INT>& sentinels) const
    {return grid.Get_MAC_Grid().Domain_Indices()+sentinels;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,grid);Initialize();}

//#####################################################################
    void Initialize();
    template<class T_ARRAYS> void Put(const T_ARRAYS& local_data,const RANGE<TV_INT>& region,T_ARRAYS& global_data) const;
    template<class T_ARRAYS> void Put(const T_ARRAYS& local_data,T_ARRAYS& global_data,const RANGE<TV_INT>& sentinels=RANGE<TV_INT>::Zero_Box()) const;
    template<class T_ARRAYS> void Get(const T_ARRAYS& global_data,T_ARRAYS& local_data) const;
    template<class T_FACE_ARRAYS> void Put_Faces(const T_FACE_ARRAYS& local_data,T_FACE_ARRAYS& global_data) const;
    template<class T_ARRAYS> T Maximum_Error(const T_ARRAYS& local_data,const T_ARRAYS& global_data,const int bandwidth,TV_INT& index,const RANGE<TV_INT>& sentinels=RANGE<TV_INT>::Zero_Box()) const;
    template<class T_FACE_ARRAYS> T Maximum_Error(const std::string& prefix,const T_FACE_ARRAYS& local_data,const T_FACE_ARRAYS& global_data,const int bandwidth,const T threshold=1e-7);
//#####################################################################
};
}
#endif
