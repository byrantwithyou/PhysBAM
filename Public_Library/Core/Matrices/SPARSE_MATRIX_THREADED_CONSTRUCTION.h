//#####################################################################
// Copyright 2017 Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_THREADED_CONSTRUCTION
//#####################################################################
#ifndef __SPARSE_MATRIX_THREADED_CONSTRUCTION__
#define __SPARSE_MATRIX_THREADED_CONSTRUCTION__

#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
namespace PhysBAM{

// User should fill in M.n; this class does not touch it.
template<class T> class SPARSE_MATRIX_THREADED_CONSTRUCTION
{
public:
    SPARSE_MATRIX_FLAT_MXN<T>& M;
    ARRAY<SPARSE_MATRIX_ENTRY<T> > A;
    ARRAY<int>& shared_num_rows;
    ARRAY<int>& shared_num_entries;
    ARRAY<int> entries_per_row;
    int tid;
    int threads;
    
    SPARSE_MATRIX_THREADED_CONSTRUCTION(SPARSE_MATRIX_FLAT_MXN<T>& M,ARRAY<int>& tmp0,ARRAY<int>& tmp1);

    void Start_Row()
    {entries_per_row.Append(0);}

    void Add_Entry(int col, T val)
    {entries_per_row.Last()++;A.Append(SPARSE_MATRIX_ENTRY<T>(col,val));}

    void Finish();

//#####################################################################
};
}
#endif
