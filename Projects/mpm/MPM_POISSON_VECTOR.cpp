//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include "MPM_POISSON_VECTOR.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_POISSON_VECTOR<TV>::
MPM_POISSON_VECTOR()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_POISSON_VECTOR<TV>::
~MPM_POISSON_VECTOR()
{
}
//#####################################################################
// Function =
//#####################################################################
template<class TV> const MPM_POISSON_VECTOR<TV>& MPM_POISSON_VECTOR<TV>::
operator= (const MPM_POISSON_VECTOR& bv)
{
    v=bv.v;
    return *this;
}
//#####################################################################
// Function +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_POISSON_VECTOR<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v+=debug_cast<const MPM_POISSON_VECTOR&>(bv).v;
    return *this;
}
//#####################################################################
// Function -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_POISSON_VECTOR<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v-=debug_cast<const MPM_POISSON_VECTOR&>(bv).v;
    return *this;
}
//#####################################################################
// Function *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_POISSON_VECTOR<TV>::
operator*=(const T a)
{
    v*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_POISSON_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    v.Copy(c,debug_cast<const MPM_POISSON_VECTOR&>(bv).v);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_POISSON_VECTOR<TV>::
Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    v.Copy(c1,debug_cast<const MPM_POISSON_VECTOR&>(bv1).v,debug_cast<const MPM_POISSON_VECTOR&>(bv2).v);
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int MPM_POISSON_VECTOR<TV>::
Raw_Size() const
{
    return Value(v.array.m)*TV::m;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& MPM_POISSON_VECTOR<TV>::
Raw_Get(int i)
{
    return v.array(i/TV::m)(i%TV::m);
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* MPM_POISSON_VECTOR<TV>::
Clone_Default() const
{
    MPM_POISSON_VECTOR<TV>* c=new MPM_POISSON_VECTOR<TV>;
    c->v.Resize(v.domain);
    return c;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void MPM_POISSON_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& w)
{
    v.Resize(debug_cast<const MPM_POISSON_VECTOR<TV>&>(w).v.domain);
}
namespace PhysBAM{
template class MPM_POISSON_VECTOR<VECTOR<float,2> >;
template class MPM_POISSON_VECTOR<VECTOR<float,3> >;
template class MPM_POISSON_VECTOR<VECTOR<double,2> >;
template class MPM_POISSON_VECTOR<VECTOR<double,3> >;
}
