//#####################################################################
// Copyright 2016, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Hybrid_Methods/Projection/MPM_PROJECTION_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PROJECTION_VECTOR<TV>::
MPM_PROJECTION_VECTOR()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PROJECTION_VECTOR<TV>::
~MPM_PROJECTION_VECTOR()
{
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_PROJECTION_VECTOR<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v+=debug_cast<const MPM_PROJECTION_VECTOR&>(bv).v;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_PROJECTION_VECTOR<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v-=debug_cast<const MPM_PROJECTION_VECTOR&>(bv).v;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& MPM_PROJECTION_VECTOR<TV>::
operator*=(const T a)
{
    v*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_PROJECTION_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    v.Copy(c,debug_cast<const MPM_PROJECTION_VECTOR&>(bv).v);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void MPM_PROJECTION_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    v.Copy(c,debug_cast<const MPM_PROJECTION_VECTOR&>(bv1).v,debug_cast<const MPM_PROJECTION_VECTOR&>(bv2).v);
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int MPM_PROJECTION_VECTOR<TV>::
Raw_Size() const
{
    return v.m;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& MPM_PROJECTION_VECTOR<TV>::
Raw_Get(int i)
{
    return v(i);
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* MPM_PROJECTION_VECTOR<TV>::
Clone_Default() const
{
    MPM_PROJECTION_VECTOR<TV>* c=new MPM_PROJECTION_VECTOR<TV>;
    c->v.Resize(v.m);
    return c;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void MPM_PROJECTION_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& vb)
{
    v.Resize(debug_cast<const MPM_PROJECTION_VECTOR<TV>&>(vb).v.m);
}
namespace PhysBAM{
template class MPM_PROJECTION_VECTOR<VECTOR<float,3> >;
template class MPM_PROJECTION_VECTOR<VECTOR<float,2> >;
template class MPM_PROJECTION_VECTOR<VECTOR<float,1> >;
template class MPM_PROJECTION_VECTOR<VECTOR<double,3> >;
template class MPM_PROJECTION_VECTOR<VECTOR<double,2> >;
template class MPM_PROJECTION_VECTOR<VECTOR<double,1> >;
}
