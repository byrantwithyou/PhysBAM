//#####################################################################
// Copyright 2009, Craig Schroeder,Russell Howes.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KRYLOV_VECTOR_CONDENSED_POISSON
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Finite_Elements/KRYLOV_VECTOR_CONDENSED_POISSON.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
KRYLOV_VECTOR_CONDENSED_POISSON()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
~KRYLOV_VECTOR_CONDENSED_POISSON()
{
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
operator+=(const KRYLOV_VECTOR_BASE<typename TV::SCALAR>& bv)
{
    v+=dynamic_cast<const KRYLOV_VECTOR_CONDENSED_POISSON&>(bv).v;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
operator-=(const KRYLOV_VECTOR_BASE<typename TV::SCALAR>& bv)
{
    v-=dynamic_cast<const KRYLOV_VECTOR_CONDENSED_POISSON&>(bv).v;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
operator*=(const T a)
{
    v*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<typename TV::SCALAR>& bv)
{
    const KRYLOV_VECTOR_CONDENSED_POISSON& v1=debug_cast<const KRYLOV_VECTOR_CONDENSED_POISSON&>(bv);
#pragma omp parallel for
        for(int i=0;i<v.m;i++)
            v(i)=c*v1.v(i);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<typename TV::SCALAR>& bv1,const KRYLOV_VECTOR_BASE<typename TV::SCALAR>& bv2)
{
    const KRYLOV_VECTOR_CONDENSED_POISSON& v1=debug_cast<const KRYLOV_VECTOR_CONDENSED_POISSON&>(bv1);
    const KRYLOV_VECTOR_CONDENSED_POISSON& v2=debug_cast<const KRYLOV_VECTOR_CONDENSED_POISSON&>(bv2);
#pragma omp parallel for
        for(int i=0;i<v.m;i++)
            v(i)=c*v1.v(i)+v2.v(i);
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
Raw_Size() const
{
    return v.m;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
Raw_Get(int i)
{
    if(i<v.m)return v(i);
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
Clone_Default() const
{
    KRYLOV_VECTOR_CONDENSED_POISSON<TV>* w=new KRYLOV_VECTOR_CONDENSED_POISSON<TV>;
    w->v.Resize(v.m);
    return w;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
Resize(const KRYLOV_VECTOR_BASE<typename TV::SCALAR>& w)
{
    v.Resize(v.m);
}
//#####################################################################
// Function Dot
//#####################################################################
template<class TV> typename TV::SCALAR KRYLOV_VECTOR_CONDENSED_POISSON<TV>::
Dot(const KRYLOV_VECTOR_BASE<typename TV::SCALAR>& bw) const
{
    const KRYLOV_VECTOR_CONDENSED_POISSON& w=debug_cast<const KRYLOV_VECTOR_CONDENSED_POISSON&>(bw);
    T result=0;
#ifdef USE_OPENMP
    int number_of_threads;
#pragma omp parallel    
#pragma omp master    
    number_of_threads=omp_get_num_threads();
    ARRAY<T> result_per_thread(number_of_threads);
    result_per_thread.Fill(0);
#pragma omp parallel for
        for(int i=0;i<v.m;i++){
            const int tid=omp_get_thread_num();
            result_per_thread(tid)+=v(i)*w.v(i);}
    for(int tid=0;tid<number_of_threads;tid++)
        result+=result_per_thread(tid);
#else
        for(int i=0;i<v.m;i++)
            result+=v(i)*w.v(i);
#endif
    return result;
}
//#####################################################################
// Function Max_Abs
//#####################################################################
template<class TV> typename TV::SCALAR KRYLOV_VECTOR_CONDENSED_POISSON<TV>::Max_Abs() const
{
#ifdef USE_OPENMP
    int number_of_threads;
#pragma omp parallel    
#pragma omp master    
    number_of_threads=omp_get_num_threads();
    ARRAY<T> result_per_thread(number_of_threads);
    result_per_thread.Fill(0);
#pragma omp parallel for
        for(int i=0;i<v.m;i++){
            const int tid=omp_get_thread_num();
            T& result=result_per_thread(tid);
            T current=fabs(v(i));
            if(current>result) result=current;}
    T max_abs=0;
    for(int tid=0;tid<number_of_threads;tid++)
        max_abs=max(max_abs,result_per_thread(tid));
#else
        T max_abs=v.Max_Abs();
#endif
    return max_abs;
}

namespace PhysBAM{
    template class KRYLOV_VECTOR_CONDENSED_POISSON<VECTOR<float,1> >;
    template class KRYLOV_VECTOR_CONDENSED_POISSON<VECTOR<float,2> >;
    template class KRYLOV_VECTOR_CONDENSED_POISSON<VECTOR<float,3> >;
    template class KRYLOV_VECTOR_CONDENSED_POISSON<VECTOR<double,1> >;
    template class KRYLOV_VECTOR_CONDENSED_POISSON<VECTOR<double,2> >;
    template class KRYLOV_VECTOR_CONDENSED_POISSON<VECTOR<double,3> >;
}
