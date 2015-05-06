//#####################################################################
// Copyright 2015, Greg Klar, Andre Pradhana
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Utilities/DEBUG_CAST.h>
#include <Hybrid_Methods/System/FLUID_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_KRYLOV_VECTOR<TV>::
FLUID_KRYLOV_VECTOR(ARRAY<int>& valid_indices)
    :valid_indices(valid_indices)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_KRYLOV_VECTOR<TV>::
~FLUID_KRYLOV_VECTOR()
{
}
//#####################################################################
// Function =
//#####################################################################
template<class TV> const FLUID_KRYLOV_VECTOR<TV>& FLUID_KRYLOV_VECTOR<TV>::
operator= (const FLUID_KRYLOV_VECTOR& bv)
{
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        p.array(i)=bv.p.array(i);}
    return *this;
}
//#####################################################################
// Function +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& FLUID_KRYLOV_VECTOR<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    const FLUID_KRYLOV_VECTOR& v=debug_cast<const FLUID_KRYLOV_VECTOR&>(bv);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        p.array(i)+=v.p.array(i);}
    return *this;
}
//#####################################################################
// Function -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& FLUID_KRYLOV_VECTOR<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    const FLUID_KRYLOV_VECTOR& v=debug_cast<const FLUID_KRYLOV_VECTOR&>(bv);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        p.array(i)-=v.p.array(i);}
    return *this;
}
//#####################################################################
// Function *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& FLUID_KRYLOV_VECTOR<TV>::
operator*=(const T a)
{
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        p.array(i)*=a;}
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void FLUID_KRYLOV_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    const FLUID_KRYLOV_VECTOR& v=debug_cast<const FLUID_KRYLOV_VECTOR&>(bv);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        p.array(i)=c*v.p.array(i);}
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void FLUID_KRYLOV_VECTOR<TV>::
Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    const FLUID_KRYLOV_VECTOR& v1=debug_cast<const FLUID_KRYLOV_VECTOR&>(bv1);
    const FLUID_KRYLOV_VECTOR& v2=debug_cast<const FLUID_KRYLOV_VECTOR&>(bv2);
#pragma omp parallel for
    for(int k=0;k<valid_indices.m;k++){
        int i=valid_indices(k);
        p.array(i)=c1*v1.p.array(i)+v2.p.array(i);}
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int FLUID_KRYLOV_VECTOR<TV>::
Raw_Size() const
{
    return valid_indices.m*TV::dimension;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& FLUID_KRYLOV_VECTOR<TV>::
Raw_Get(int i)
{
    return p.array(valid_indices(i));
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* FLUID_KRYLOV_VECTOR<TV>::
Clone_Default() const
{
    FLUID_KRYLOV_VECTOR<TV>* c=new FLUID_KRYLOV_VECTOR<TV>(valid_indices);
    c->p.Resize(p.domain);
    return c;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void FLUID_KRYLOV_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& w)
{
    p.Resize(debug_cast<const FLUID_KRYLOV_VECTOR<TV>&>(w).p.domain);
}
namespace PhysBAM{
template class FLUID_KRYLOV_VECTOR<VECTOR<float,2> >;
template class FLUID_KRYLOV_VECTOR<VECTOR<float,3> >;
template class FLUID_KRYLOV_VECTOR<VECTOR<double,2> >;
template class FLUID_KRYLOV_VECTOR<VECTOR<double,3> >;
}
