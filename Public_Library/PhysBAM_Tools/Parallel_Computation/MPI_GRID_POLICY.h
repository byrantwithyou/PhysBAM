//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPI_GRID_POLICY__
#define __MPI_GRID_POLICY__

namespace PhysBAM{

template<class T_GRID> class MPI_UNIFORM_GRID;

template<class T_GRID>
struct MPI_GRID_POLICY
{
    typedef MPI_UNIFORM_GRID<T_GRID> MPI_GRID;
    typedef T_GRID PARALLEL_GRID;
};
}
#endif
