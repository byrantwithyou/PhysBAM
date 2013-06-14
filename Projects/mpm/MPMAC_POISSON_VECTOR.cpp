//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include "MPMAC_POISSON_VECTOR.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPMAC_POISSON_VECTOR<TV>::
MPMAC_POISSON_VECTOR()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPMAC_POISSON_VECTOR<TV>::
~MPMAC_POISSON_VECTOR()
{
}
//#####################################################################
// Function =
//#####################################################################
template<class TV> const MPMAC_POISSON_VECTOR<TV>& MPMAC_POISSON_VECTOR<TV>::
operator= (const MPMAC_POISSON_VECTOR& bv)
{
    v=bv.v;
    return *this;
}
//#####################################################################
// Function +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPMAC_POISSON_VECTOR<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v+=debug_cast<const MPMAC_POISSON_VECTOR&>(bv).v;
    return *this;
}
//#####################################################################
// Function -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPMAC_POISSON_VECTOR<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v-=debug_cast<const MPMAC_POISSON_VECTOR&>(bv).v;
    return *this;
}
//#####################################################################
// Function *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPMAC_POISSON_VECTOR<TV>::
operator*=(const T a)
{
    v*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPMAC_POISSON_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    v.Copy(c,debug_cast<const MPMAC_POISSON_VECTOR&>(bv).v);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPMAC_POISSON_VECTOR<TV>::
Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    v.Copy(c1,debug_cast<const MPMAC_POISSON_VECTOR&>(bv1).v,debug_cast<const MPMAC_POISSON_VECTOR&>(bv2).v);
}
//#####################################################################
// Function Dot
//#####################################################################
template<class TV> typename TV::SCALAR MPMAC_POISSON_VECTOR<TV>::
Dot(const KRYLOV_VECTOR_BASE<T>& bv) const
{
    return v.array.Dot(debug_cast<const MPMAC_POISSON_VECTOR&>(bv).v.array);
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int MPMAC_POISSON_VECTOR<TV>::
Raw_Size() const
{
    return Value(v.array.m);
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& MPMAC_POISSON_VECTOR<TV>::
Raw_Get(int i)
{
    return v.array(i);
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* MPMAC_POISSON_VECTOR<TV>::
Clone_Default() const
{
    MPMAC_POISSON_VECTOR<TV>* c=new MPMAC_POISSON_VECTOR<TV>;
    c->v.Resize(v.domain);
    return c;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void MPMAC_POISSON_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& w)
{
    v.Resize(debug_cast<const MPMAC_POISSON_VECTOR<TV>&>(w).v.domain);
}
namespace PhysBAM{
template class MPMAC_POISSON_VECTOR<VECTOR<float,2> >;
template class MPMAC_POISSON_VECTOR<VECTOR<float,3> >;
template class MPMAC_POISSON_VECTOR<VECTOR<double,2> >;
template class MPMAC_POISSON_VECTOR<VECTOR<double,3> >;
}
