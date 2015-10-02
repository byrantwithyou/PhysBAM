//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar, Chenfanfu Jiang, Chuyuan Fu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Utilities/DEBUG_CAST.h>
#include <Hybrid_Methods/System/MPM_KKT_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_KKT_KRYLOV_VECTOR<TV>::
MPM_KKT_KRYLOV_VECTOR(ARRAY<int>& valid_indices,ARRAY<int>& valid_p_indices)
    :valid_indices(valid_indices),valid_p_indices(valid_p_indices)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_KKT_KRYLOV_VECTOR<TV>::
~MPM_KKT_KRYLOV_VECTOR()
{
}
//#####################################################################
// Function =
//#####################################################################
template<class TV> const MPM_KKT_KRYLOV_VECTOR<TV>& MPM_KKT_KRYLOV_VECTOR<TV>::
operator= (const MPM_KKT_KRYLOV_VECTOR& bv)
{
    PHYSBAM_ASSERT(bv.lambda.m==lambda.m);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)=bv.u.array(i);}
#pragma omp parallel for
    for(int k=0;k<valid_p_indices.m;k++){
        int i=valid_p_indices(k);
        p.array(i)=bv.p.array(i);}
    for (int k=0;k<lambda.m;k++)
        lambda(k)=bv.lambda(k);
    return *this;
}
//#####################################################################
// Function +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_KKT_KRYLOV_VECTOR<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    const MPM_KKT_KRYLOV_VECTOR& v=debug_cast<const MPM_KKT_KRYLOV_VECTOR&>(bv);
    PHYSBAM_ASSERT(v.lambda.m==lambda.m);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)+=v.u.array(i);}
#pragma omp parallel for
    for(int k=0;k<valid_p_indices.m;k++){
        int i=valid_p_indices(k);
        p.array(i)+=v.p.array(i);}
    for (int k=0;k<lambda.m;k++)
        lambda(k)+=v.lambda(k);
    return *this;
}
//#####################################################################
// Function -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_KKT_KRYLOV_VECTOR<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    const MPM_KKT_KRYLOV_VECTOR& v=debug_cast<const MPM_KKT_KRYLOV_VECTOR&>(bv);
    PHYSBAM_ASSERT(v.lambda.m==lambda.m);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)-=v.u.array(i);}
#pragma omp parallel for
    for(int k=0;k<valid_p_indices.m;k++){
        int i=valid_p_indices(k);
        p.array(i)-=v.p.array(i);}
    for (int k=0;k<lambda.m;k++)
        lambda(k)-=v.lambda(k);
    return *this;
}
//#####################################################################
// Function *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_KKT_KRYLOV_VECTOR<TV>::
operator*=(const T a)
{
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)*=a;}
#pragma omp parallel for
    for(int k=0;k<valid_p_indices.m;k++){
        int i=valid_p_indices(k);
        p.array(i)*=a;}
    for (int k=0;k<lambda.m;k++)
        lambda(k)*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_KKT_KRYLOV_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    const MPM_KKT_KRYLOV_VECTOR& v=debug_cast<const MPM_KKT_KRYLOV_VECTOR&>(bv);
    PHYSBAM_ASSERT(v.lambda.m==lambda.m);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)=c*v.u.array(i);}
#pragma omp parallel for
    for(int k=0;k<valid_p_indices.m;k++){
        int i=valid_p_indices(k);
        p.array(i)=c*v.p.array(i);}
    for (int k=0;k<lambda.m;k++)
        lambda(k)=c*v.lambda(k);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_KKT_KRYLOV_VECTOR<TV>::
Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    const MPM_KKT_KRYLOV_VECTOR& v1=debug_cast<const MPM_KKT_KRYLOV_VECTOR&>(bv1);
    const MPM_KKT_KRYLOV_VECTOR& v2=debug_cast<const MPM_KKT_KRYLOV_VECTOR&>(bv2);
    PHYSBAM_ASSERT(v1.lambda.m==v2.lambda.m);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)=c1*v1.u.array(i)+v2.u.array(i);}
    for(int k=0;k<valid_p_indices.m;k++){
        int i=valid_p_indices(k);
        p.array(i)=c1*v1.p.array(i)+v2.p.array(i);}
    for (int k=0;k<lambda.m;k++)
        lambda(k)=c1*v1.lambda(k)+v2.lambda(k);
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int MPM_KKT_KRYLOV_VECTOR<TV>::
Raw_Size() const
{
    return (TV::m)*valid_indices.m+valid_p_indices.m+lambda.m;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& MPM_KKT_KRYLOV_VECTOR<TV>::
Raw_Get(int i)
{
    int velocity_size=TV::m*valid_indices.m;
    if(i<velocity_size) return u.array(valid_indices(i/TV::m))(i%TV::m);
    else if(i<velocity_size+valid_p_indices.m) return p.array(valid_p_indices(i-velocity_size));
    else return lambda(i-velocity_size-valid_p_indices.m);
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* MPM_KKT_KRYLOV_VECTOR<TV>::
Clone_Default() const
{
    MPM_KKT_KRYLOV_VECTOR<TV>* c=new MPM_KKT_KRYLOV_VECTOR<TV>(valid_indices,valid_p_indices);
    c->u.Resize(u.domain);
    c->p.Resize(p.domain);
    c->lambda.Resize(lambda.m);
    return c;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void MPM_KKT_KRYLOV_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& w)
{
    u.Resize(debug_cast<const MPM_KKT_KRYLOV_VECTOR<TV>&>(w).u.domain);
    p.Resize(debug_cast<const MPM_KKT_KRYLOV_VECTOR<TV>&>(w).p.domain);
    lambda.Resize(debug_cast<const MPM_KKT_KRYLOV_VECTOR<TV>&>(w).lambda.m);
}
namespace PhysBAM{
template class MPM_KKT_KRYLOV_VECTOR<VECTOR<float,1> >;
template class MPM_KKT_KRYLOV_VECTOR<VECTOR<float,2> >;
template class MPM_KKT_KRYLOV_VECTOR<VECTOR<float,3> >;
template class MPM_KKT_KRYLOV_VECTOR<VECTOR<double,1> >;
template class MPM_KKT_KRYLOV_VECTOR<VECTOR<double,2> >;
template class MPM_KKT_KRYLOV_VECTOR<VECTOR<double,3> >;
}
