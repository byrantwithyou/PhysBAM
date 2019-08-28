//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __KRYLOV_SOLVER_FEM__
#define __KRYLOV_SOLVER_FEM__
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#if USE_LAPACK
#include <cblas.h>
#include <lapacke.h>
#elif USE_MKL
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#include <mkl_spblas.h>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace PhysBAM{

template<typename T>
struct KRYLOV_VECTOR_FEM:public KRYLOV_VECTOR_BASE<T>
{
    typedef ARRAY<T> TV;
    TV v;
    bool deep_copy;
    int threads=1;

    KRYLOV_VECTOR_FEM(int th):deep_copy(false),threads(th)
    {
    }
    KRYLOV_VECTOR_FEM(TV vector,int th):v(vector),deep_copy(false),threads(th)
    {
    }
    virtual ~KRYLOV_VECTOR_FEM()
    {
    }
    static void Scalar_Times_Vector(double a,const ARRAY<double>& x,ARRAY<double>& y,int threads)
    {
#if 0
        int s=x.m/threads,l=x.m-s*(threads-1);
#pragma omp parallel for
        for(int tid=0;tid<threads;tid++)
        {
            int n=tid==threads-1?l:s;
            cblas_daxpy(n,a,x.base_pointer+tid*s,1,y.base_pointer+tid*s,1);
        }
#else
        cblas_daxpy(x.m,a,x.base_pointer,1,y.base_pointer,1);
#endif
    }
    static void Scalar_Times_Vector(float a,const ARRAY<float>& x,ARRAY<float>& y,int threads)
    {
#if 0
        int s=x.m/threads,l=x.m-s*(threads-1);
#pragma omp parallel for
        for(int tid=0;tid<threads;tid++)
        {
            int n=tid==threads-1?l:s;
            cblas_saxpy(n,a,x.base_pointer+tid*s,1,y.base_pointer+tid*s,1);
        }
#else
        cblas_saxpy(x.m,a,x.base_pointer,1,y.base_pointer,1);
#endif
    }
    static void Scale_Vector(double a,ARRAY<double>& x,int threads)
    {
#if 0
        int s=x.m/threads,l=x.m-s*(threads-1);
#pragma omp parallel for
        for(int tid=0;tid<threads;tid++)
        {
            int n=tid==threads-1?l:s;
            cblas_dscal(n,a,x.base_pointer+tid*s,1);
        }
#else
        cblas_dscal(x.m,a,x.base_pointer,1);
#endif
    }
    static void Scale_Vector(float a,ARRAY<float >& x,int threads)
    {
#if 0
        int s=x.m/threads,l=x.m-s*(threads-1);
#pragma omp parallel for
        for(int tid=0;tid<threads;tid++)
        {
            int n=tid==threads-1?l:s;
            cblas_sscal(n,a,x.base_pointer+tid*s,1);
        }
#else
        cblas_sscal(x.m,a,x.base_pointer,1);
#endif
    }
    KRYLOV_VECTOR_BASE<T>& operator+=(const KRYLOV_VECTOR_BASE<T>& bv) override
    {
        const ARRAY<T>& w=dynamic_cast<const KRYLOV_VECTOR_FEM<T>&>(bv).v;
        Scalar_Times_Vector((T)1,w,v,threads);
        return *this;
    }
    KRYLOV_VECTOR_BASE<T>& operator-=(const KRYLOV_VECTOR_BASE<T>& bv) override
    {
        const ARRAY<T>& w=dynamic_cast<const KRYLOV_VECTOR_FEM<T>&>(bv).v;
        Scalar_Times_Vector(-(T)1,w,v,threads);
        return *this;
    }
    KRYLOV_VECTOR_BASE<T>& operator*=(const T a) override
    {
        Scale_Vector(a,v,threads);
        return *this;
    }
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv) override
    {
        const ARRAY<T>& w=dynamic_cast<const KRYLOV_VECTOR_FEM<T>&>(bv).v;
        int s=v.m/threads,l=v.m-s*(threads-1);
#pragma omp parallel for
        for(int tid=0;tid<threads;tid++)
        {
            int n=tid==threads-1?l:s;
            for(int i=0,j=tid*s;i<n;i++,j++)
                v(j)=c*w(j);
        }
    }
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) override
    {
        const ARRAY<T>& w1=dynamic_cast<const KRYLOV_VECTOR_FEM<T>&>(bv1).v;
        const ARRAY<T>& w2=dynamic_cast<const KRYLOV_VECTOR_FEM<T>&>(bv2).v;
        int s=v.m/threads,l=v.m-s*(threads-1);
#pragma omp parallel for
        for(int tid=0;tid<threads;tid++)
        {
            int n=tid==threads-1?l:s;
            for(int i=0,j=tid*s;i<n;i++,j++)
                v(j)=c*w1(j)+w2(j);
        }
    }
    int Raw_Size() const override
    {
        return v.Size();
    }
    T& Raw_Get(int i) override
    {
        return v(i);
    }
    KRYLOV_VECTOR_BASE<T>* Clone_Default() const override
    {
        KRYLOV_VECTOR_FEM<T>* c=new KRYLOV_VECTOR_FEM<T>(threads);
        c->v.Resize(v.Size());
        c->deep_copy=true;
        return c;
    }
    void Resize(const KRYLOV_VECTOR_BASE<T>& u) override
    {
        if(!deep_copy) return;
        const KRYLOV_VECTOR_FEM<T>& w=debug_cast<const KRYLOV_VECTOR_FEM<T>&>(u);
        v.Resize(w.v.Size());
    }
    void Get(ARRAY_VIEW<T> a) const override
    {
        a=v;
    }
    void Set(ARRAY_VIEW<const T> a) override
    {
        v=a;
    }
};

template<typename T>
struct KRYLOV_SYSTEM_FEM:public KRYLOV_SYSTEM_BASE<T>
{
    typedef KRYLOV_VECTOR_FEM<T> VECTOR_T;
    typedef SPARSE_MATRIX_FLAT_MXN<T> T_MATRIX;
    typedef T_MATRIX T_MATRIX_PRECON;
    int threads=1;
    const T_MATRIX& A;
    const T_MATRIX_PRECON* P;
    mutable VECTOR_T* temp_vector;
    ARRAY<VECTOR_T*> nullspace_vectors; // must be orthogonal and normalized
    bool use_mkl_sparse=false;
#if USE_MKL
    ARRAY<int> ja;
    ARRAY<T> entries;
    sparse_matrix_t spA;
    matrix_descr descA;
#endif

    KRYLOV_SYSTEM_FEM(const T_MATRIX& A_input,int th,bool use_mkl_sp):
        KRYLOV_SYSTEM_BASE<T>(false,true),threads(th),A(A_input),use_mkl_sparse(use_mkl_sp)
    {
#ifdef USE_OPENMP
        omp_set_num_threads(threads);
#endif

#if USE_MKL
        if(use_mkl_sparse)
        {
            mkl_set_num_threads(threads);
            ja.Resize(A.A.m);
            entries.Resize(A.A.m);
            for(int i=0;i<A.A.m;i++)
            {
                ja(i)=A.A(i).j;
                entries(i)=A.A(i).a;
            }
            Init_Sparse_Matrix((T)0);
            descA.type=SPARSE_MATRIX_TYPE_GENERAL;
            descA.diag=SPARSE_DIAG_NON_UNIT;
            const int EST_ITERATIONS=100000;
            mkl_sparse_set_mv_hint(spA,SPARSE_OPERATION_NON_TRANSPOSE,descA,EST_ITERATIONS);
            mkl_sparse_optimize(spA);
        }
#endif
    }
    virtual ~KRYLOV_SYSTEM_FEM() = default;

    void Set_Preconditioner(const T_MATRIX_PRECON& preconditioner,VECTOR_T& vector)
    {
        this->use_preconditioner=true;
        P=&preconditioner;
        temp_vector=&vector;
    }
#if USE_MKL
    void Init_Sparse_Matrix(double dummy)
    {
        mkl_sparse_d_create_csr(&spA,SPARSE_INDEX_BASE_ZERO,A.m,A.n,
            A.offsets.base_pointer,A.offsets.base_pointer+1,
            ja.base_pointer,entries.base_pointer);
    }
    void Init_Sparse_Matrix(float dummy)
    {
        mkl_sparse_s_create_csr(&spA,SPARSE_INDEX_BASE_ZERO,A.m,A.n,
            A.offsets.base_pointer,A.offsets.base_pointer+1,
            ja.base_pointer,entries.base_pointer);
    }
    void Matrix_Times_Vector(const ARRAY<double>& x,ARRAY<double>& result) const
    {
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,1,
            spA,descA,(T*)x.base_pointer,0,result.base_pointer);
    }
    void Matrix_Times_Vector(const ARRAY<float>& x,ARRAY<float>& result) const
    {
        mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE,1,
            spA,descA,(T*)x.base_pointer,0,result.base_pointer);
    }
#endif
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result,bool transpose=false) const override
    {
        const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x);
        VECTOR_T& vresult=dynamic_cast<VECTOR_T&>(result);
        if(use_mkl_sparse)
        {
#if USE_MKL
            vresult.v.Resize(vx.v.m);
            Matrix_Times_Vector(vx.v,vresult.v);
#else
            PHYSBAM_FATAL_ERROR();
#endif
        }
        else
        {
            A.Times_Threaded(vx.v,vresult.v);
        }
    }
    static float Dot_Product(const ARRAY<float>& u,const ARRAY<float>& v,int threads)
    {
        PHYSBAM_ASSERT(u.m==v.m);
#if 0
        int s=u.m/threads,l=u.m-s*(threads-1);
        float r=0;
#pragma omp parallel for reduction(+:r)
        for(int tid=0;tid<threads;tid++)
        {
            float d=cblas_sdot(tid==threads-1?l:s,u.base_pointer+tid*s,1,v.base_pointer+tid*s,1);
            r=r+d;
        }
        return r;
#else
        return cblas_sdot(u.m,u.base_pointer,1,v.base_pointer,1);
#endif
    }
    static double Dot_Product(const ARRAY<double>& u,const ARRAY<double>& v,int threads)
    {
        PHYSBAM_ASSERT(u.m==v.m);
#if 0
        int s=u.m/threads,l=u.m-s*(threads-1);
        double r=0;
#pragma omp parallel for reduction(+:r)
        for(int tid=0;tid<threads;tid++)
        {
            double d=cblas_ddot(tid==threads-1?l:s,u.base_pointer+tid*s,1,v.base_pointer+tid*s,1);
            r=r+d;
        }
        return r;
#else
        return cblas_ddot(u.m,u.base_pointer,1,v.base_pointer,1);
#endif
    }
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const override
    {
        const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x),vy=dynamic_cast<const VECTOR_T&>(y);
        return Dot_Product(vx.v,vy.v,threads);
    }
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const override
    {
        const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x);
        int s=vx.v.m/threads,l=vx.v.m-s*(threads-1);
        T r=0;
#pragma omp parallel for reduction(max:r)
        for(int tid=0;tid<threads;tid++)
        {
            double d=0;
            int n=tid==threads-1?l:s;
            for(int i=0,j=tid*s;i<n;i++,j++)
            {
                d=max(d,abs(vx.v(j)));
            }
            r=max(r,d);
        }
        return r;
    }
    void Project(KRYLOV_VECTOR_BASE<T>& x) const override
    {
        for(auto p:nullspace_vectors)
            x.Copy(-Inner_Product(x,*p),*p,x);
    }
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const override {}
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const override
    {
        Project(x);
    } 
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const override
    {
        const VECTOR_T& vr=dynamic_cast<const VECTOR_T&>(r);
        VECTOR_T& vz=dynamic_cast<VECTOR_T&>(z);
        P->Solve_Forward_Substitution(vr.v,temp_vector->v,true);
        P->Solve_Backward_Substitution(temp_vector->v,vz.v,false,true);
    }
};
}
#endif
