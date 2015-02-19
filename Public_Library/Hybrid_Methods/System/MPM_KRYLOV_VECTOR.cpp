//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Utilities/DEBUG_CAST.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_KRYLOV_VECTOR<TV>::
MPM_KRYLOV_VECTOR(ARRAY<int>& valid_indices)
    :valid_indices(valid_indices)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_KRYLOV_VECTOR<TV>::
~MPM_KRYLOV_VECTOR()
{
}
//#####################################################################
// Function =
//#####################################################################
template<class TV> const MPM_KRYLOV_VECTOR<TV>& MPM_KRYLOV_VECTOR<TV>::
operator= (const MPM_KRYLOV_VECTOR& bv)
{
    u.array.Subset(valid_indices)=bv.u.array.Subset(valid_indices);
    return *this;
}
//#####################################################################
// Function +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_KRYLOV_VECTOR<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    u.array.Subset(valid_indices)+=debug_cast<const MPM_KRYLOV_VECTOR&>(bv).u.array.Subset(valid_indices);
    return *this;
}
//#####################################################################
// Function -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_KRYLOV_VECTOR<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    u.array.Subset(valid_indices)-=debug_cast<const MPM_KRYLOV_VECTOR&>(bv).u.array.Subset(valid_indices);
    return *this;
}
//#####################################################################
// Function *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_KRYLOV_VECTOR<TV>::
operator*=(const T a)
{
    u.array.Subset(valid_indices)*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_KRYLOV_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    u.array.Subset(valid_indices)=c*debug_cast<const MPM_KRYLOV_VECTOR&>(bv).u.array.Subset(valid_indices);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_KRYLOV_VECTOR<TV>::
Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    u.array.Subset(valid_indices)=c1*debug_cast<const MPM_KRYLOV_VECTOR&>(bv1).u.array.Subset(valid_indices)+
        debug_cast<const MPM_KRYLOV_VECTOR&>(bv2).u.array.Subset(valid_indices);
}
//#####################################################################
// Function Dot
//#####################################################################
template<class TV> typename TV::SCALAR MPM_KRYLOV_VECTOR<TV>::
Dot(const KRYLOV_VECTOR_BASE<T>& bv) const
{
    return u.array.Subset(valid_indices).Dot(debug_cast<const MPM_KRYLOV_VECTOR&>(bv).u.array.Subset(valid_indices));
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int MPM_KRYLOV_VECTOR<TV>::
Raw_Size() const
{
    return valid_indices.m*TV::dimension;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& MPM_KRYLOV_VECTOR<TV>::
Raw_Get(int i)
{
    return u.array(valid_indices(i/TV::m))(i%TV::m);
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* MPM_KRYLOV_VECTOR<TV>::
Clone_Default() const
{
    MPM_KRYLOV_VECTOR<TV>* c=new MPM_KRYLOV_VECTOR<TV>(valid_indices);
    c->u.Resize(u.domain);
    return c;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void MPM_KRYLOV_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& w)
{
    u.Resize(debug_cast<const MPM_KRYLOV_VECTOR<TV>&>(w).u.domain);
}
namespace PhysBAM{
template class MPM_KRYLOV_VECTOR<VECTOR<float,2> >;
template class MPM_KRYLOV_VECTOR<VECTOR<float,3> >;
template class MPM_KRYLOV_VECTOR<VECTOR<double,2> >;
template class MPM_KRYLOV_VECTOR<VECTOR<double,3> >;
}
