//#####################################################################
// Copyright 2017 Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class APPEND_HOLDER
//#####################################################################
#ifndef __APPEND_HOLDER__
#define __APPEND_HOLDER__

#include <Core/Arrays/ARRAY.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{

template<class T>
class APPEND_HOLDER
{
public:
    struct SPACER
    {
        ARRAY<T> array;
        int offset;
        char spacer[128-sizeof(ARRAY<T>)-sizeof(int)];
    };
    ARRAY<T>& orig_array;
    ARRAY<SPACER> thread_array;
    int threads;

    APPEND_HOLDER(ARRAY<T>& orig_array)
        :orig_array(orig_array)
    {
    }

    // Run from within #pragma omp parallel
    void Init()
    {
#ifdef USE_OPENMP
        threads=omp_get_num_threads();
        if(threads>1)
            thread_array.Resize(threads);
#else
        threads=1;
#endif
    }

    // Run from within #pragma omp parallel
    ARRAY<T>& Array()
    {
#ifdef USE_OPENMP
        if(!thread_array.m) return orig_array;
        return thread_array(omp_get_thread_num()).array;
#else
        return orig_array;
#endif
    }

    void Combine()
    {
        if(!thread_array.m) return;
        int size=0;
        for(int i=0;i<threads;i++){
            thread_array(i).offset=size;
            size+=thread_array(i).array.m;}
        orig_array.Resize(size,false,false);
#pragma omp parallel
        {
#ifdef USE_OPENMP
            int t=omp_get_thread_num();
#else
            int t=0;
#endif
            const ARRAY<T>& p=thread_array(t).array;
            for(int i=0,a=p.m,b=thread_array(t).offset;i<a;i++)
                orig_array(i+b)=p(i);
        }
    }
};
}
#endif
