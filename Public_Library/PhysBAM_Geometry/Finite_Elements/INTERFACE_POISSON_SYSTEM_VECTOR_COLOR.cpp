//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_VECTOR_COLOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
INTERFACE_POISSON_SYSTEM_VECTOR_COLOR()
{
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
    u=v.u;
    q=v.q;
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
operator+=(const BASE& bv)
{
    const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR& v=debug_cast<const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR&>(bv);
    u+=v.u;
    q+=v.q;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
operator-=(const BASE& bv)
{
    const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR& v=debug_cast<const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR&>(bv);
    u-=v.u;
    q-=v.q;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
operator*=(const T a)
{
    for(int c=0;c<colors;c++)
        u(c)*=a;
    q*=a;
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
        u(c).Copy(c1,v1.u(c));
    q.Copy(c1,v1.q);
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
        u(c).Copy(c1,v1.u(c),v2.u(c));
    q.Copy(c1,v1.q,v2.q);
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Print() const
{
    // Flat print
    for(int c=0;c<colors;c++)
        for(int k=0;k<u(c).n;k++)
            LOG::cout<<u(c)(k)<<" ";
    for(int k=0;k<q.n;k++)
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
        size+=u(c).n;
    size+=q.n;
    return size;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Raw_Get(int i)
{
    for(int c=0;c<colors;c++){
        if(i<u(c).n) return u(c)(i);
            i-=u(c).n;}
    if(i<q.n) return q(i);
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
    for(int c=0;c<colors;c++)
        v->u(c).Resize(u(c).n);
    v->q.Resize(q.n);
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
    for(int c=0;c<colors;c++)
        u(c).Resize(cs.u(c).n);
    q.Resize(cs.q.n);
}
//#####################################################################
// Function Dot
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Dot(const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>& v) const
{
    T dot=0;
    for(int c=0;c<colors;c++)
        dot+=u(c).Dot(v.u(c));
    dot+=q.Dot(v.q);
    return dot;
}
//#####################################################################
// Function Magnitude_Squared
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Magnitude_Squared() const
{
    return Dot(*this);
}
//#####################################################################
// Function Magnitude
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Magnitude() const
{
    return sqrt(this->Magnitude_Squared());
}
//#####################################################################
// Function Max_Abs
//#####################################################################
template<class TV> typename TV::SCALAR INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Max_Abs() const
{
    T max_abs=0;
    for(int c=0;c<colors;c++)
        max_abs=max(u(c).Max_Abs(),max_abs);
    max_abs=max(q.Max_Abs(),max_abs);
    return max_abs;
}
//#####################################################################
// Function Normalize
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Normalize()
{
    (*this)*=1/this->Magnitude();
}
//#####################################################################
// Function Scale
//#####################################################################
template<class TV> void INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>::
Scale(const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>& v,const INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<TV>& s)
{
    for(int c=0;c<colors;c++)
        for(int k=0;k<u(c).n;k++)
            u(c)(k)=v.u(c)(k)*s.u(c)(k);
    for(int k=0;k<q.n;k++)
        q(k)=v.q(k)*s.q(k);
}
//#####################################################################
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<float,1> >;
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<float,2> >;
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<double,1> >;
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<double,2> >;
template class INTERFACE_POISSON_SYSTEM_VECTOR_COLOR<VECTOR<double,3> >;
#endif
