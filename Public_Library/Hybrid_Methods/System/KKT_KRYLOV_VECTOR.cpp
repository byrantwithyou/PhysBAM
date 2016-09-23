//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar, Chenfanfu Jiang, Chuyuan Fu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Utilities/DEBUG_CAST.h>
#include <Hybrid_Methods/System/KKT_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> KKT_KRYLOV_VECTOR<TV>::
KKT_KRYLOV_VECTOR(ARRAY<int>& valid_indices)
    :valid_indices(valid_indices)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> KKT_KRYLOV_VECTOR<TV>::
~KKT_KRYLOV_VECTOR()
{
}
//#####################################################################
// Function =
//#####################################################################
template<class TV> const KKT_KRYLOV_VECTOR<TV>& KKT_KRYLOV_VECTOR<TV>::
operator= (const KKT_KRYLOV_VECTOR& bv)
{
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)=bv.u.array(i);
        p.array(i)=bv.p.array(i);}
    return *this;
}
//#####################################################################
// Function +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& KKT_KRYLOV_VECTOR<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    const KKT_KRYLOV_VECTOR& v=debug_cast<const KKT_KRYLOV_VECTOR&>(bv);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)+=v.u.array(i);
        p.array(i)+=v.p.array(i);}
    return *this;
}
//#####################################################################
// Function -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& KKT_KRYLOV_VECTOR<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    const KKT_KRYLOV_VECTOR& v=debug_cast<const KKT_KRYLOV_VECTOR&>(bv);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)-=v.u.array(i);
        p.array(i)-=v.p.array(i);}
    return *this;
}
//#####################################################################
// Function *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& KKT_KRYLOV_VECTOR<TV>::
operator*=(const T a)
{
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)*=a;
        p.array(i)*=a;}
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void KKT_KRYLOV_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    const KKT_KRYLOV_VECTOR& v=debug_cast<const KKT_KRYLOV_VECTOR&>(bv);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)=c*v.u.array(i);
        p.array(i)=c*v.p.array(i);}
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void KKT_KRYLOV_VECTOR<TV>::
Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    const KKT_KRYLOV_VECTOR& v1=debug_cast<const KKT_KRYLOV_VECTOR&>(bv1);
    const KKT_KRYLOV_VECTOR& v2=debug_cast<const KKT_KRYLOV_VECTOR&>(bv2);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)=c1*v1.u.array(i)+v2.u.array(i);
        p.array(i)=c1*v1.p.array(i)+v2.p.array(i);}
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int KKT_KRYLOV_VECTOR<TV>::
Raw_Size() const
{
    return (TV::m+1)*valid_indices.m;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& KKT_KRYLOV_VECTOR<TV>::
Raw_Get(int i)
{
    int velocity_size=TV::m*valid_indices.m;
    if(i<velocity_size) return u.array(valid_indices(i/TV::m))(i%TV::m);
    else return p.array(valid_indices(i-velocity_size));
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* KKT_KRYLOV_VECTOR<TV>::
Clone_Default() const
{
    KKT_KRYLOV_VECTOR<TV>* c=new KKT_KRYLOV_VECTOR<TV>(valid_indices);
    c->u.Resize(u.domain);
    c->p.Resize(p.domain);
    return c;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void KKT_KRYLOV_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& w)
{
    u.Resize(debug_cast<const KKT_KRYLOV_VECTOR<TV>&>(w).u.domain);
    p.Resize(debug_cast<const KKT_KRYLOV_VECTOR<TV>&>(w).p.domain);
}
namespace PhysBAM{
template class KKT_KRYLOV_VECTOR<VECTOR<float,2> >;
template class KKT_KRYLOV_VECTOR<VECTOR<float,3> >;
template class KKT_KRYLOV_VECTOR<VECTOR<double,2> >;
template class KKT_KRYLOV_VECTOR<VECTOR<double,3> >;
}
