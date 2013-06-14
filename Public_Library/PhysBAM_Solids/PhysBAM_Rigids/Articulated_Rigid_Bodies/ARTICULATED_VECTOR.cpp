//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> ARTICULATED_VECTOR<TV>::
ARTICULATED_VECTOR()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ARTICULATED_VECTOR<TV>::
~ARTICULATED_VECTOR()
{
}
//#####################################################################
// Function =
//#####################################################################
template<class TV> const ARTICULATED_VECTOR<TV>& ARTICULATED_VECTOR<TV>::
operator= (const ARTICULATED_VECTOR& bv)
{
    v=bv.v;
    return *this;
}
//#####################################################################
// Function +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& ARTICULATED_VECTOR<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v+=debug_cast<const ARTICULATED_VECTOR&>(bv).v;
    return *this;
}
//#####################################################################
// Function -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& ARTICULATED_VECTOR<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v-=debug_cast<const ARTICULATED_VECTOR&>(bv).v;
    return *this;
}
//#####################################################################
// Function *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& ARTICULATED_VECTOR<TV>::
operator*=(const T a)
{
    v*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void ARTICULATED_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    v=c*debug_cast<const ARTICULATED_VECTOR&>(bv).v;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void ARTICULATED_VECTOR<TV>::
Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    v=c1*debug_cast<const ARTICULATED_VECTOR&>(bv1).v+debug_cast<const ARTICULATED_VECTOR&>(bv2).v;
}
//#####################################################################
// Function Dot
//#####################################################################
template<class TV> typename TV::SCALAR ARTICULATED_VECTOR<TV>::
Dot(const KRYLOV_VECTOR_BASE<T>& bv) const
{
    return v.Dot(debug_cast<const ARTICULATED_VECTOR&>(bv).v);
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int ARTICULATED_VECTOR<TV>::
Raw_Size() const
{
    return Value(v.m)*TWIST<TV>::dimension;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& ARTICULATED_VECTOR<TV>::
Raw_Get(int i)
{
    int o=i%TWIST<TV>::dimension,n=i/TWIST<TV>::dimension;
    if(o<TV::dimension) return v(JOINT_ID(n)).linear(o);
    return v(JOINT_ID(n)).angular(o-TV::dimension);
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* ARTICULATED_VECTOR<TV>::
Clone_Default() const
{
    ARTICULATED_VECTOR<TV>* c=new ARTICULATED_VECTOR<TV>;
    c->v.Resize(v.m);
    return c;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void ARTICULATED_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& w)
{
    v.Resize(debug_cast<const ARTICULATED_VECTOR<TV>&>(w).v.m);
}
namespace PhysBAM{
template class ARTICULATED_VECTOR<VECTOR<float,1> >;
template class ARTICULATED_VECTOR<VECTOR<float,2> >;
template class ARTICULATED_VECTOR<VECTOR<float,3> >;
template class ARTICULATED_VECTOR<VECTOR<double,1> >;
template class ARTICULATED_VECTOR<VECTOR<double,2> >;
template class ARTICULATED_VECTOR<VECTOR<double,3> >;
}
