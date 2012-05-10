//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM_VECTOR_NEW
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR_NEW.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
INTERFACE_STOKES_SYSTEM_VECTOR_NEW()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
~INTERFACE_STOKES_SYSTEM_VECTOR_NEW()
{
}
//#####################################################################
// Operator =
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
operator=(const INTERFACE_STOKES_SYSTEM_VECTOR_NEW& v)
{
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s]=v.u(i)[s];
    for(int s=0;s<2;s++) p[s]=v.p[s];
    q=v.q;
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
operator+=(const BASE& bv)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_NEW& v=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_NEW&>(bv);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s]+=v.u(i)[s];
    for(int s=0;s<2;s++) p[s]+=v.p[s];
    q+=v.q;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
operator-=(const BASE& bv)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_NEW& v=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_NEW&>(bv);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s]-=v.u(i)[s];
    for(int s=0;s<2;s++) p[s]-=v.p[s];
    q-=v.q;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
operator*=(const T a)
{
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s]*=a;
    for(int s=0;s<2;s++) p[s]*=a;;
    q*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Copy(const T c,const BASE& bv)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_NEW& v=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_NEW&>(bv);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            VECTOR_ND<T>::Copy(c,v.u(i)[s],u(i)[s]);
    for(int s=0;s<2;s++) VECTOR_ND<T>::Copy(c,v.p[s],p[s]);
    VECTOR_ND<T>::Copy(c,v.q,q);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_NEW& v1=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_NEW&>(bv1);
    const INTERFACE_STOKES_SYSTEM_VECTOR_NEW& v2=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_NEW&>(bv2);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            VECTOR_ND<T>::Copy(c1,v1.u(i)[s],v2.u(i)[s],u(i)[s]);
    for(int s=0;s<2;s++) VECTOR_ND<T>::Copy(c1,v1.p[s],v2.p[s],p[s]);
    VECTOR_ND<T>::Copy(c1,v1.q,v2.q,q);
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Print() const
{
    // Flat print
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            for(int k=0;k<u(i)[s].n;k++)
                LOG::cout<<u(i)[s](k)<<" ";
    for(int s=0;s<2;s++)
        for(int k=0;k<p[s].n;k++)
            LOG::cout<<p[s](k)<<" ";
    for(int k=0;k<q.n;k++)
        LOG::cout<<q(k)<<" ";
    LOG::cout<<std::endl;
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Raw_Size() const
{
    int size=0;
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            size+=u(i)[s].n;
    for(int s=0;s<2;s++) size+=p[s].n;
    size+=q.n;
    return size;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Raw_Get(int i)
{
    for(int k=0;k<TV::m;k++)
        for(int s=0;s<2;s++){
            if(i<u(k)[s].n) return u(k)[s](i);
            i-=u(k)[s].n;}
    for(int s=0;s<2;s++){
        if(i<p[s].n) return p[s](i);
        i-=p[s].n;}
    if(i<q.n) return q(i);
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Clone_Default() const
{
    INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>* v=new INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>;
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            v->u(i)[s].Resize(u(i)[s].n);
    for(int s=0;s<2;s++) v->p[s].Resize(p[s].n);
    v->q.Resize(q.n);
    return v;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& v)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& cs=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>&>(v);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s].Resize(cs.u(i)[s].n);
    for(int s=0;s<2;s++) p[s].Resize(cs.p[s].n);
    q.Resize(cs.q.n);
}
//#####################################################################
// Function Dot
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Dot(const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& v) const
{
    T dot=0;
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            dot+=u(i)[s].Dot(v.u(i)[s]);
    for(int s=0;s<2;s++) dot+=p[s].Dot(v.p[s]);
    dot+=q.Dot(v.q);
    return dot;
}
//#####################################################################
// Function Magnitude_Squared
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Magnitude_Squared() const
{
    T ms=0;
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            ms+=u(i)[s].Magnitude_Squared();
    for(int s=0;s<2;s++) ms+=p[s].Magnitude_Squared();
    ms+=q.Magnitude_Squared();
    return ms;
}
//#####################################################################
// Function Magnitude
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Magnitude() const
{
    return sqrt(this->Magnitude_Squared());
}
//#####################################################################
// Function Max_Abs
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Max_Abs() const
{
    T max_abs=0;
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            max_abs=max(u(i)[s].Max_Abs(),max_abs);
    for(int s=0;s<2;s++) max_abs=max(p[s].Max_Abs(),max_abs);
    max_abs=max(q.Max_Abs(),max_abs);
    return max_abs;
}
//#####################################################################
// Function Normalize
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Normalize()
{
    (*this)*=1/this->Magnitude();
}
//#####################################################################
// Function Scale
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Scale(const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& v,const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& c)
{
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            for(int k=0;k<u(i)[s].n;k++)
                u(i)[s](k)=v.u(i)[s](k)*c.u(i)[s](k);
    for(int s=0;s<2;s++)
        for(int k=0;k<p[s].n;k++)
            p[s](k)=v.p[s](k)*c.p[s](k);
    for(int k=0;k<q.n;k++)
        q(k)=v.q(k)*c.q(k);
}
//#####################################################################
// Function Scale
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>::
Scale(const INTERFACE_STOKES_SYSTEM_VECTOR_NEW<TV>& c)
{
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            for(int k=0;k<u(i)[s].n;k++)
                u(i)[s](k)*=c.u(i)[s](k);
    for(int s=0;s<2;s++)
        for(int k=0;k<p[s].n;k++)
            p[s](k)*=c.p[s](k);
    for(int k=0;k<q.n;k++)
        q(k)*=c.q(k);
}
//#####################################################################
template class INTERFACE_STOKES_SYSTEM_VECTOR_NEW<VECTOR<float,1> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR_NEW<VECTOR<float,2> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR_NEW<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INTERFACE_STOKES_SYSTEM_VECTOR_NEW<VECTOR<double,1> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR_NEW<VECTOR<double,2> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR_NEW<VECTOR<double,3> >;
#endif
