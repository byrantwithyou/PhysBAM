//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM_VECTOR_COLOR
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR_COLOR.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
INTERFACE_STOKES_SYSTEM_VECTOR_COLOR()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
~INTERFACE_STOKES_SYSTEM_VECTOR_COLOR()
{
}
//#####################################################################
// Operator =
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
operator=(const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR& v)
{
#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<colors;c++)
#pragma omp task
                u(i)(c)=v.u(i)(c);
        for(int c=0;c<colors;c++)
#pragma omp task
            p(c)=v.p(c);
#pragma omp task
        q=v.q;
    }
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
operator+=(const BASE& bv)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR& v=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR&>(bv);
#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<colors;c++)
#pragma omp task
                u(i)(c)+=v.u(i)(c);
        for(int c=0;c<colors;c++)
#pragma omp task
            p(c)+=v.p(c);
#pragma omp task
        q+=v.q;
    }
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
operator-=(const BASE& bv)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR& v=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR&>(bv);
#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<colors;c++)
#pragma omp task
                u(i)(c)-=v.u(i)(c);
        for(int c=0;c<colors;c++)
#pragma omp task
            p(c)-=v.p(c);
#pragma omp task
        q-=v.q;
    }
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
operator*=(const T a)
{
#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<colors;c++)
#pragma omp task
                u(i)(c)*=a;
        for(int c=0;c<colors;c++)
#pragma omp task
            p(c)*=a;;
#pragma omp task
        q*=a;
    }
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Copy(const T c1,const BASE& bv1)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR& v1=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR&>(bv1);
#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<colors;c++)
#pragma omp task
                u(i)(c).Copy(c1,v1.u(i)(c));
        for(int c=0;c<colors;c++)
#pragma omp task
            p(c).Copy(c1,v1.p(c));
#pragma omp task
        q.Copy(c1,v1.q);
    }
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR& v1=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR&>(bv1);
    const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR& v2=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR&>(bv2);
#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<colors;c++)
#pragma omp task
                u(i)(c).Copy(c1,v1.u(i)(c),v2.u(i)(c));
        for(int c=0;c<colors;c++)
#pragma omp task
            p(c).Copy(c1,v1.p(c),v2.p(c));
#pragma omp task
        q.Copy(c1,v1.q,v2.q);
    }
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Print() const
{
    // Flat print
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<colors;c++)
            for(int k=0;k<u(i)(c).m;k++)
                LOG::cout<<u(i)(c)(k)<<" ";
    for(int c=0;c<colors;c++)
        for(int k=0;k<p(c).m;k++)
            LOG::cout<<p(c)(k)<<" ";
    for(int k=0;k<q.m;k++)
        LOG::cout<<q(k)<<" ";
    LOG::cout<<std::endl;
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Raw_Size() const
{
    int size=0;
    for(int i=0;i<TV::m;i++)
        for(int c=0;c<colors;c++)
            size+=u(i)(c).m;
    for(int c=0;c<colors;c++) size+=p(c).m;
    size+=q.m;
    return size;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Raw_Get(int i)
{
    for(int k=0;k<TV::m;k++)
        for(int c=0;c<colors;c++){
            if(i<u(k)(c).m) return u(k)(c)(i);
            i-=u(k)(c).m;}
    for(int c=0;c<colors;c++){
        if(i<p(c).m) return p(c)(i);
        i-=p(c).m;}
    if(i<q.m) return q(i);
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Clone_Default() const
{
    INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>* v=new INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>;
    v->colors=colors;
    for(int i=0;i<TV::m;i++){
        v->u(i).Resize(colors);
        for(int c=0;c<colors;c++)
            v->u(i)(c).Resize(u(i)(c).m);}
    v->p.Resize(colors);
    for(int c=0;c<colors;c++) v->p(c).Resize(p(c).m);
    v->q.Resize(q.m);
    return v;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& v)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& cs=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>&>(v);
    colors=cs.colors;
    for(int i=0;i<TV::m;i++){
        u(i).Resize(colors);
        for(int c=0;c<colors;c++)
            u(i)(c).Resize(cs.u(i)(c).m);}
    p.Resize(colors);
    for(int c=0;c<colors;c++) p(c).Resize(cs.p(c).m);
    q.Resize(cs.q.m);
}
//#####################################################################
// Function Dot
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Dot(const KRYLOV_VECTOR_BASE<T>& bv) const
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& v=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>&>(bv);
    T dot=0;
#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<colors;c++)
#pragma omp task
            {
                const T tmp=u(i)(c).Dot(v.u(i)(c));
#pragma omp critical
                dot+=tmp;
            }
        for(int c=0;c<colors;c++)
#pragma omp task
        {
            const T tmp=p(c).Dot(v.p(c));
#pragma omp critical
            dot+=tmp;
        }
#pragma omp task
        {
            const T tmp=q.Dot(v.q);
#pragma omp critical
            dot+=tmp;
        }
    }
    return dot;
}
//#####################################################################
// Function Magnitude_Squared
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Magnitude_Squared() const
{
    return Dot(*this);
}
//#####################################################################
// Function Magnitude
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Magnitude() const
{
    return sqrt(this->Magnitude_Squared());
}
//#####################################################################
// Function Max_Abs
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Max_Abs() const
{
    T max_abs=0;
#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<colors;c++)
#pragma omp task
            {
                const T tmp=u(i)(c).Max_Abs();
#pragma omp critical
                max_abs=max(tmp,max_abs);
            }
        for(int c=0;c<colors;c++)
#pragma omp task
        {
            const T tmp=p(c).Max_Abs();
#pragma omp critical
            max_abs=max(tmp,max_abs);
        }
#pragma omp task
        {
            const T tmp=q.Max_Abs();
#pragma omp critical
            max_abs=max(tmp,max_abs);
        }
    }
    return max_abs;
}
//#####################################################################
// Function Normalize
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Normalize()
{
    (*this)*=1/this->Magnitude();
}
//#####################################################################
// Function Scale
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>::
Scale(const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& v,const INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV>& s)
{
#pragma omp parallel
#pragma omp single
    {
        for(int i=0;i<TV::m;i++)
            for(int c=0;c<colors;c++)
#pragma omp task
                for(int k=0;k<u(i)(c).m;k++)
                    u(i)(c)(k)=v.u(i)(c)(k)*s.u(i)(c)(k);
        for(int c=0;c<colors;c++)
#pragma omp task
            for(int k=0;k<p(c).m;k++)
                p(c)(k)=v.p(c)(k)*s.p(c)(k);
#pragma omp task
        for(int k=0;k<q.m;k++)
            q(k)=v.q(k)*s.q(k);
    }
}
//#####################################################################
namespace PhysBAM{
template class INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<VECTOR<float,1> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<VECTOR<float,2> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<VECTOR<float,3> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<VECTOR<double,1> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<VECTOR<double,2> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<VECTOR<double,3> >;
}
