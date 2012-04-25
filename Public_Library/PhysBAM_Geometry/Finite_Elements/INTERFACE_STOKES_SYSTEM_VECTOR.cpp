//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM_VECTOR
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
INTERFACE_STOKES_SYSTEM_VECTOR()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
~INTERFACE_STOKES_SYSTEM_VECTOR()
{
}
//#####################################################################
// Operator =
//#####################################################################
template<class TV> INTERFACE_STOKES_SYSTEM_VECTOR<TV>& INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
operator=(const INTERFACE_STOKES_SYSTEM_VECTOR& v)
{
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s]=v.u(i)[s];
    for(int s=0;s<2;s++) p[s]=v.p[s];
    for(int i=0;i<TV::m;i++) q(i)=v.q(i);
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
operator+=(const BASE& bv)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR& v=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR&>(bv);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s]+=v.u(i)[s];
    for(int s=0;s<2;s++) p[s]+=v.p[s];
    for(int i=0;i<TV::m;i++) q(i)+=v.q(i);
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
operator-=(const BASE& bv)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR& v=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR&>(bv);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s]-=v.u(i)[s];
    for(int s=0;s<2;s++) p[s]-=v.p[s];
    for(int i=0;i<TV::m;i++) q(i)-=v.q(i);
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
operator*=(const T a)
{
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s]*=a;
    for(int s=0;s<2;s++) p[s]*=a;;
    for(int i=0;i<TV::m;i++) q(i)*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
Copy(const T c,const BASE& bv)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR& v=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR&>(bv);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            VECTOR_ND<T>::Copy(c,v.u(i)[s],u(i)[s]);
    for(int s=0;s<2;s++) VECTOR_ND<T>::Copy(c,v.p[s],p[s]);
    for(int i=0;i<TV::m;i++) VECTOR_ND<T>::Copy(c,v.q(i),q(i));
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR& v1=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR&>(bv1);
    const INTERFACE_STOKES_SYSTEM_VECTOR& v2=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR&>(bv2);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            VECTOR_ND<T>::Copy(c1,v1.u(i)[s],v2.u(i)[s],u(i)[s]);
    for(int s=0;s<2;s++) VECTOR_ND<T>::Copy(c1,v1.p[s],v2.p[s],p[s]);
    for(int i=0;i<TV::m;i++) VECTOR_ND<T>::Copy(c1,v1.q(i),v2.q(i),q(i));
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
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
    for(int i=0;i<TV::m;i++)
        for(int k=0;k<q(i).n;k++)
            LOG::cout<<q(i)(k)<<" ";
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
Raw_Size() const
{
    int size=0;
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            size+=u(i)[s].n;
    for(int s=0;s<2;s++) size+=p[s].n;
    for(int i=0;i<TV::m;i++) size+=q(i).n;
    return size;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
Raw_Get(int i)
{
    for(int k=0;k<TV::m;k++)
        for(int s=0;s<2;s++){
            if(i<u(k)[s].n) return u(k)[s](i);
            i-=u(k)[s].n;}
    for(int s=0;s<2;s++){
        if(i<p[s].n) return p[s](i);
        i-=p[s].n;}
    for(int k=0;k<TV::m;k++){
        if(i<q(k).n) return q(k)(i);
        i-=q(k).n;}
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
Clone_Default() const
{
    INTERFACE_STOKES_SYSTEM_VECTOR<TV>* v=new INTERFACE_STOKES_SYSTEM_VECTOR<TV>;
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            v->u(i)[s].Resize(u(i)[s].n);
    for(int s=0;s<2;s++) v->p[s].Resize(p[s].n);
    for(int i=0;i<TV::m;i++) v->q(i).Resize(q(i).n);
    return v;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& v)
{
    const INTERFACE_STOKES_SYSTEM_VECTOR<TV>& cs=debug_cast<const INTERFACE_STOKES_SYSTEM_VECTOR<TV>&>(v);
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            u(i)[s].Resize(cs.u(i)[s].n);
    for(int s=0;s<2;s++) p[s].Resize(cs.p[s].n);
    for(int i=0;i<TV::m;i++) q(i).Resize(cs.q(i).n);

}
//#####################################################################
template class INTERFACE_STOKES_SYSTEM_VECTOR<VECTOR<float,1> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR<VECTOR<float,2> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INTERFACE_STOKES_SYSTEM_VECTOR<VECTOR<double,1> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR<VECTOR<double,2> >;
template class INTERFACE_STOKES_SYSTEM_VECTOR<VECTOR<double,3> >;
#endif
