//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_SYSTEM_VECTOR
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
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
        for(int s=0;s<2;s++){
            u[i][s]=v.u[i][s];
            q[i][s]=v.q[i][s];}
    for(int s=0;s<2;s++) p[s]=v.p[s];
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
        for(int s=0;s<2;s++){
            u[i][s]+=v.u[i][s];
            q[i][s]+=v.q[i][s];}
    for(int s=0;s<2;s++) p[s]+=v.p[s];
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
        for(int s=0;s<2;s++){
            u[i][s]-=v.u[i][s];
            q[i][s]-=v.q[i][s];}
    for(int s=0;s<2;s++) p[s]-=v.p[s];
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
operator*=(const T a)
{
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            u[i][s]*=a;
            q[i][s]*=a;}
    for(int s=0;s<2;s++) p[s]*=a;;
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
        for(int s=0;s<2;s++){
            ARRAY<T>::Copy(c,v.u[i][s],u[i][s]);
            ARRAY<T>::Copy(c,v.q[i][s],q[i][s]);}
    for(int s=0;s<2;s++) ARRAY<T>::Copy(c,v.p[s],p[s]);
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
        for(int s=0;s<2;s++){
            ARRAY<T>::Copy(c1,v1.u[i][s],v2.u[i][s],u[i][s]);
            ARRAY<T>::Copy(c1,v1.q[i][s],v2.q[i][s],q[i][s]);}
    for(int s=0;s<2;s++) ARRAY<T>::Copy(c1,v1.p[s],v2.p[s],p[s]);
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
            for(int k=0;k<u[i][s].m;k++)
                LOG::cout<<u[i][s](k)<<" ";
    for(int s=0;s<2;s++)
        for(int k=0;k<p[s].m;k++)
            LOG::cout<<p[s](k)<<" ";
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++)
            for(int k=0;k<q[i][s].m;k++)
                LOG::cout<<q[i][s](k)<<" ";
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
            size+=u[i][s].m+q[i][s].m;
    for(int s=0;s<2;s++) size+=p[s].m;
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
            if(i<u[k][s].m) return u[k][s](i);
            i-=u[k][s].m;}
    for(int s=0;s<2;s++){
        if(i<p[s].m) return p[s](i);
        i-=p[s].m;}
    for(int k=0;k<TV::m;k++)
        for(int s=0;s<2;s++){
            if(i<q[k][s].m) return q[k][s](i);
            i-=q[k][s].m;}
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* INTERFACE_STOKES_SYSTEM_VECTOR<TV>::
Clone_Default() const
{
    INTERFACE_STOKES_SYSTEM_VECTOR<TV>* v=new INTERFACE_STOKES_SYSTEM_VECTOR<TV>;
    for(int i=0;i<TV::m;i++)
        for(int s=0;s<2;s++){
            v->u[i][s].Resize(u[i][s].m);
            v->q[i][s].Resize(q[i][s].m);}
    for(int s=0;s<2;s++) v->p[s].Resize(p[s].m);
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
        for(int s=0;s<2;s++){
            u[i][s].Resize(cs.u[i][s].m);
            q[i][s].Resize(cs.q[i][s].m);}
    for(int s=0;s<2;s++) p[s].Resize(cs.p[s].m);
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
