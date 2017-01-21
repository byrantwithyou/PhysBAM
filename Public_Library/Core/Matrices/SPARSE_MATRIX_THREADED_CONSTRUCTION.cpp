//#####################################################################
// Copyright 2017 Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/SPARSE_MATRIX_THREADED_CONSTRUCTION.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
// User should fill in M.n; this class does not touch it.
template<class T> SPARSE_MATRIX_THREADED_CONSTRUCTION<T>::
SPARSE_MATRIX_THREADED_CONSTRUCTION(SPARSE_MATRIX_FLAT_MXN<T>& M,ARRAY<int>& tmp0,ARRAY<int>& tmp1)
    :M(M),shared_num_rows(tmp0),shared_num_entries(tmp1),tid(0),threads(1)
{
#ifdef USE_OPENMP
    threads=omp_get_num_threads();
    tid=omp_get_thread_num();
#endif
#pragma omp single
    {
      shared_num_rows.Resize(threads+1);
      shared_num_entries.Resize(threads+1);
    }
#pragma omp barrier
}
//#####################################################################
// Function Start_Second_Pass
//#####################################################################
template<class T> void SPARSE_MATRIX_THREADED_CONSTRUCTION<T>::
Finish()
{
    shared_num_rows(tid)=entries_per_row.m;
    shared_num_entries(tid)=entries_per_row.Sum();
#pragma omp barrier
#pragma omp single
    {
        int last_row=0,last_entries=0;
        for(int i=0;i<threads;i++){
            int s=shared_num_rows(i);
            shared_num_rows(i)=last_row;
            last_row+=s;
            int e=shared_num_entries(i);
            shared_num_entries(i)=last_entries;
            last_entries+=e;}
        shared_num_rows.Last()=last_row;
        shared_num_entries.Last()=last_entries;
        M.m=last_row;
        M.offsets.Resize(last_row+1,false);
        M.A.Resize(last_entries,false);
    }
#pragma omp barrier
    int first_row=shared_num_rows(tid),last=shared_num_entries(tid);
    for(int i=0;i<A.m;i++)
        M.A(i+last)=A(i);
    for(int i=0;i<entries_per_row.m;i++){
        M.offsets(i+first_row)=last;
        M.A.Array_View(last,entries_per_row(i)).Sort();
        last+=entries_per_row(i);}
    if(tid==threads-1) M.offsets.Last()=last;
#pragma omp barrier
}
template class SPARSE_MATRIX_THREADED_CONSTRUCTION<float>;
template class SPARSE_MATRIX_THREADED_CONSTRUCTION<double>;
}
