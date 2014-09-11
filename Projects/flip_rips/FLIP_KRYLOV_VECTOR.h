//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLIP_KRYLOV_VECTOR__
#define __FLIP_KRYLOV_VECTOR__

#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#ifdef USE_OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num(){return 0;}
#endif

namespace PhysBAM{

template<class TV>
class FLIP_KRYLOV_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef KRYLOV_VECTOR_BASE<T> BASE;

public:

    ARRAY<T,TV_INT> p;
    const ARRAY<TV_INT>& interior_cells;
    const int threads;

    FLIP_KRYLOV_VECTOR(const ARRAY<TV_INT>& interior_cells_input,const int threads)
        :interior_cells(interior_cells_input),threads(threads){}
    virtual ~FLIP_KRYLOV_VECTOR(){}

    FLIP_KRYLOV_VECTOR& operator=(const FLIP_KRYLOV_VECTOR& v)
    {
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            p(index)=v.p(index);}
        return *this;
    }

    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE
    {
        const FLIP_KRYLOV_VECTOR<TV>& v=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(bv);
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            p(index)+=v.p(index);}
        return *this;
    }

    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE
    {
        const FLIP_KRYLOV_VECTOR<TV>& v=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(bv);
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            p(index)-=v.p(index);}
        return *this;
    }

    BASE& operator*=(const T a) PHYSBAM_OVERRIDE
    {
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            p(index)*=a;}
        return *this;
    }

    void Copy(const T c,const BASE& bv) PHYSBAM_OVERRIDE
    {
        const FLIP_KRYLOV_VECTOR<TV>& v=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(bv);
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            p(index)=c*v.p(index);}
    }

    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE
    {
        const FLIP_KRYLOV_VECTOR<TV>& v1=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(bv1);
        const FLIP_KRYLOV_VECTOR<TV>& v2=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(bv2);
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            p(index)=c1*v1.p(index)+v2.p(index);}
    }

    T Dot(const KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE
    {
        const FLIP_KRYLOV_VECTOR<TV>& v=debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(bv);

        ARRAY<double> result_per_thread(threads);
#pragma omp parallel for
        for(int k=0;k<interior_cells.m;k++){
            const TV_INT index=interior_cells(k);
            const int tid=omp_get_thread_num();
            result_per_thread(tid)+=(double)p(index)*(double)v.p(index);}
        double result=0;
        for(int tid=0;tid<threads;tid++)
            result+=result_per_thread(tid);
        return result;
    }

    BASE* Clone_Default() const
    {
        FLIP_KRYLOV_VECTOR<TV>* c=new FLIP_KRYLOV_VECTOR<TV>(interior_cells,threads);
        c->p.Resize(p.domain);
        return c;
    }

    void Resize(const KRYLOV_VECTOR_BASE<T>& w)
    {p.Resize(debug_cast<const FLIP_KRYLOV_VECTOR<TV>&>(w).p.domain);}

    int Raw_Size() const PHYSBAM_OVERRIDE
    {return interior_cells.m;}

    T& Raw_Get(int j) PHYSBAM_OVERRIDE
    {return p(interior_cells(j));}
};
}
#endif
