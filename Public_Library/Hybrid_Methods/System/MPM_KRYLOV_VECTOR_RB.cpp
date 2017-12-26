//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Utilities/DEBUG_CAST.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR_RB.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_KRYLOV_VECTOR_RB<TV>::
MPM_KRYLOV_VECTOR_RB(ARRAY<int>& valid_indices)
    :valid_indices(valid_indices)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_KRYLOV_VECTOR_RB<TV>::
~MPM_KRYLOV_VECTOR_RB()
{
}
//#####################################################################
// Function =
//#####################################################################
template<class TV> const MPM_KRYLOV_VECTOR_RB<TV>& MPM_KRYLOV_VECTOR_RB<TV>::
operator= (const MPM_KRYLOV_VECTOR_RB& bv)
{
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)=bv.u.array(i);}
#pragma omp parallel for
    for(int k=0;k<twists.m;k++)
        twists(k)=bv.twists(k);
    return *this;
}
//#####################################################################
// Function +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_KRYLOV_VECTOR_RB<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    const MPM_KRYLOV_VECTOR_RB& v=debug_cast<const MPM_KRYLOV_VECTOR_RB&>(bv);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)+=v.u.array(i);}
#pragma omp parallel for
    for(int k=0;k<twists.m;k++)
        twists(k)+=v.twists(k);
    return *this;
}
//#####################################################################
// Function -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_KRYLOV_VECTOR_RB<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    const MPM_KRYLOV_VECTOR_RB& v=debug_cast<const MPM_KRYLOV_VECTOR_RB&>(bv);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)-=v.u.array(i);}
#pragma omp parallel for
    for(int k=0;k<twists.m;k++)
        twists(k)-=v.twists(k);
    return *this;
}
//#####################################################################
// Function *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_KRYLOV_VECTOR_RB<TV>::
operator*=(const T a)
{
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)*=a;}
#pragma omp parallel for
    for(int k=0;k<twists.m;k++)
        twists(k)*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_KRYLOV_VECTOR_RB<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    const MPM_KRYLOV_VECTOR_RB& v=debug_cast<const MPM_KRYLOV_VECTOR_RB&>(bv);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)=c*v.u.array(i);}
#pragma omp parallel for
    for(int k=0;k<twists.m;k++)
        twists(k)=c*v.twists(k);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_KRYLOV_VECTOR_RB<TV>::
Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    const MPM_KRYLOV_VECTOR_RB& v1=debug_cast<const MPM_KRYLOV_VECTOR_RB&>(bv1);
    const MPM_KRYLOV_VECTOR_RB& v2=debug_cast<const MPM_KRYLOV_VECTOR_RB&>(bv2);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        u.array(i)=c1*v1.u.array(i)+v2.u.array(i);}
#pragma omp parallel for
    for(int k=0;k<twists.m;k++)
        twists(k)=c1*v1.twists(k)+v2.twists(k);
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int MPM_KRYLOV_VECTOR_RB<TV>::
Raw_Size() const
{
    return valid_indices.m*TV::m+twists.m*TWIST<TV>::m;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& MPM_KRYLOV_VECTOR_RB<TV>::
Raw_Get(int i)
{
    if(i<valid_indices.m*TV::m) return u.array(valid_indices(i/TV::m))(i%TV::m);
    i-=valid_indices.m*TV::m;
    TWIST<TV>& t=twists(i/TWIST<TV>::m);
    if(i<TV::m) return t.linear(i);
    return t.angular(i-TV::m);
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* MPM_KRYLOV_VECTOR_RB<TV>::
Clone_Default() const
{
    MPM_KRYLOV_VECTOR_RB<TV>* c=new MPM_KRYLOV_VECTOR_RB<TV>(valid_indices);
    c->u.Resize(u.domain);
    c->twists.Resize(twists.m);
    return c;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void MPM_KRYLOV_VECTOR_RB<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& w)
{
    const MPM_KRYLOV_VECTOR_RB<TV>& v=debug_cast<const MPM_KRYLOV_VECTOR_RB<TV>&>(w);
    u.Resize(v.u.domain);
    twists.Resize(v.twists.m);
}
namespace PhysBAM{
template class MPM_KRYLOV_VECTOR_RB<VECTOR<float,2> >;
template class MPM_KRYLOV_VECTOR_RB<VECTOR<float,3> >;
template class MPM_KRYLOV_VECTOR_RB<VECTOR<double,2> >;
template class MPM_KRYLOV_VECTOR_RB<VECTOR<double,3> >;
}
