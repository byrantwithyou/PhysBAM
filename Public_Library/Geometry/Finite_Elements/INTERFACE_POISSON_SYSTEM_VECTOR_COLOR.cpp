//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR
//#####################################################################
#include <Core/Log/LOG.h>
#include <Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_VECTOR_COLOR.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
INTERFACE_POISSON_SYSTEM_VECTOR_COLOR()
{
#ifdef USE_OPENMP
#pragma omp parallel
#pragma omp master
    threads=omp_get_num_threads();
    result_per_thread.Resize(threads);
#endif
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
~INTERFACE_POISSON_SYSTEM_VECTOR_COLOR()
{
}
//#####################################################################
// Operator =
//#####################################################################
template<class TV> INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>& INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
operator=(const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR& v)
{
    for(int c=0;c<colors;c++)
#pragma omp parallel for
        for(int i=0;i<u(c).m;i++)
            u(c)(i)=v.u(c)(i);
#pragma omp parallel for
        for(int i=0;i<q.m;i++)
            q(i)=v.q(i);
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
operator+=(const BASE& bv)
{
    const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR& v=debug_cast<const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR&>(bv);
    for(int c=0;c<colors;c++)
#pragma omp parallel for
        for(int i=0;i<u(c).m;i++)
            u(c)(i)+=v.u(c)(i);
#pragma omp parallel for
        for(int i=0;i<q.m;i++)
            q(i)+=v.q(i);
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
operator-=(const BASE& bv)
{
    const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR& v=debug_cast<const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR&>(bv);
    for(int c=0;c<colors;c++)
#pragma omp parallel for
        for(int i=0;i<u(c).m;i++)
            u(c)(i)-=v.u(c)(i);
#pragma omp parallel for
        for(int i=0;i<q.m;i++)
            q(i)-=v.q(i);
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
operator*=(const T a)
{
    for(int c=0;c<colors;c++)
#pragma omp parallel for
        for(int i=0;i<u(c).m;i++)
            u(c)(i)*=a;
#pragma omp parallel for
        for(int i=0;i<q.m;i++)
            q(i)*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Copy(const T c1,const BASE& bv1)
{
    const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR& v1=debug_cast<const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR&>(bv1);
    for(int c=0;c<colors;c++)
#pragma omp parallel for
        for(int i=0;i<u(c).m;i++)
            u(c)(i)=c1*v1.u(c)(i);
#pragma omp parallel for
        for(int i=0;i<q.m;i++)
            q(i)=c1*v1.q(i);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR& v1=debug_cast<const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR&>(bv1);
    const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR& v2=debug_cast<const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR&>(bv2);
    for(int c=0;c<colors;c++)
#pragma omp parallel for
        for(int i=0;i<u(c).m;i++)
            u(c)(i)=c1*v1.u(c)(i)+v2.u(c)(i);
#pragma omp parallel for
        for(int i=0;i<q.m;i++)
            q(i)=c1*v1.q(i)+v2.q(i);
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Print() const
{
    // Flat print
    for(int c=0;c<colors;c++)
        for(int k=0;k<u(c).m;k++)
            LOG::cout<<u(c)(k)<<" ";
    for(int k=0;k<q.m;k++)
        LOG::cout<<q(k)<<" ";
    LOG::cout<<std::endl;
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Raw_Size() const
{
    int size=0;
    for(int c=0;c<colors;c++)
        size+=u(c).m;
    size+=q.m;
    return size;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Raw_Get(int i)
{
    for(int c=0;c<colors;c++){
        if(i<u(c).m) return u(c)(i);
            i-=u(c).m;}
    if(i<q.m) return q(i);
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Clone_Default() const
{
    INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>* v=new INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>;
    v->colors=colors;
    v->u.Resize(colors);
#pragma omp parallel for
    for(int c=0;c<colors;c++)
        v->u(c).Resize(u(c).m);
    v->q.Resize(q.m);
    return v;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& v)
{
    const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>& cs=debug_cast<const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>&>(v);
    colors=cs.colors;
    u.Resize(colors);
#pragma omp parallel for
    for(int c=0;c<colors;c++)
        u(c).Resize(cs.u(c).m);
    q.Resize(cs.q.m);
}
//#####################################################################
// Function Set_Zero
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Set_Zero()
{
    for(int c=0;c<colors;c++)
#pragma omp parallel for
        for(int i=0;i<u(c).m;i++)
            u(c)(i)=0;
#pragma omp parallel for
        for(int i=0;i<q.m;i++)
            q(i)=0;
}
//#####################################################################
// Function Max_Abs
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Max_Abs() const
{
#ifdef USE_OPENMP
    result_per_thread.Fill(0);
    for(int c=0;c<colors;c++)
#pragma omp parallel for
        for(int i=0;i<u(c).m;i++){
            const int tid=omp_get_thread_num();
            T& result=result_per_thread(tid);
            T current=fabs(u(c)(i));
            if(current>result) result=current;}
#pragma omp parallel for
    for(int i=0;i<q.m;i++){
        const int tid=omp_get_thread_num();
        T& result=result_per_thread(tid);
        T current=fabs(q(i));
        if(current>result) result=current;}
    T max_abs=0;
    for(int tid=0;tid<threads;tid++)
        max_abs=max(max_abs,result_per_thread(tid));
#else
    T max_abs=q.Max_Abs();
    for(int c=0;c<colors;c++)
        max_abs=max(max_abs,u(c).Max_Abs());
#endif
    return max_abs;
}
//#####################################################################
// Function Scale
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Scale(const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>& v,const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>& s)
{
    for(int c=0;c<colors;c++)
#pragma omp parallel for
        for(int k=0;k<u(c).m;k++)
            u(c)(k)=v.u(c)(k)*s.u(c)(k);
#pragma omp parallel for
    for(int k=0;k<q.m;k++)
        q(k)=v.q(k)*s.q(k);
}
//#####################################################################
// Function Get
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Get(ARRAY_VIEW<T> a) const
{
    int k=0;
    for(const auto&i:u)
        for(const auto&j:i)
            a(k++)=j;
    for(const auto&i:q)
        a(k++)=i;
}
//#####################################################################
// Function Set
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Set(ARRAY_VIEW<const T> a)
{
    int k=0;
    for(auto&i:u)
        for(auto&j:i)
            j=a(k++);
    for(auto&i:q)
        i=a(k++);
}
//#####################################################################
namespace PhysBAM{
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<float,1> >;
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<float,2> >;
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<float,3> >;
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<double,1> >;
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<double,2> >;
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<double,3> >;
}
