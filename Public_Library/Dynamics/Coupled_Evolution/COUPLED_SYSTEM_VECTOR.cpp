//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COUPLED_SYSTEM_VECTOR
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Dynamics/Coupled_Evolution/COUPLED_SYSTEM_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COUPLED_SYSTEM_VECTOR<TV>::
COUPLED_SYSTEM_VECTOR()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COUPLED_SYSTEM_VECTOR<TV>::
~COUPLED_SYSTEM_VECTOR()
{
}
//#####################################################################
// Operator =
//#####################################################################
template<class TV> COUPLED_SYSTEM_VECTOR<TV>& COUPLED_SYSTEM_VECTOR<TV>::
operator=(const COUPLED_SYSTEM_VECTOR& v)
{
    pressure=v.pressure;
    lambda=v.lambda;
    force_coefficients=v.force_coefficients;
    viscous_force_coefficients=v.viscous_force_coefficients;
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& COUPLED_SYSTEM_VECTOR<TV>::
operator+=(const BASE& bv)
{
    const COUPLED_SYSTEM_VECTOR& v=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv);
    pressure+=v.pressure;
    lambda+=v.lambda;
    force_coefficients+=v.force_coefficients;
    viscous_force_coefficients+=v.viscous_force_coefficients;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& COUPLED_SYSTEM_VECTOR<TV>::
operator-=(const BASE& bv)
{
    const COUPLED_SYSTEM_VECTOR& v=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv);
    pressure-=v.pressure;
    lambda-=v.lambda;
    force_coefficients-=v.force_coefficients;
    viscous_force_coefficients-=v.viscous_force_coefficients;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& COUPLED_SYSTEM_VECTOR<TV>::
operator*=(const T a)
{
    pressure*=a;
    lambda*=a;
    force_coefficients*=a;
    viscous_force_coefficients*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void COUPLED_SYSTEM_VECTOR<TV>::
Copy(const T c,const BASE& bv)
{
    const COUPLED_SYSTEM_VECTOR& v=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv);
    assert(v.pressure.m==pressure.m);
    pressure.Copy(c,v.pressure);
    lambda=c*v.lambda;
    force_coefficients=c*v.force_coefficients;
    viscous_force_coefficients=c*v.viscous_force_coefficients;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void COUPLED_SYSTEM_VECTOR<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const COUPLED_SYSTEM_VECTOR& v1=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv1);
    const COUPLED_SYSTEM_VECTOR& v2=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv2);
    assert(v1.pressure.m==v2.pressure.m && pressure.m==v1.pressure.m);
    pressure.Copy(c1,v1.pressure,v2.pressure);
    lambda=c1*v1.lambda+v2.lambda;
    force_coefficients=c1*v1.force_coefficients+v2.force_coefficients;
    viscous_force_coefficients=c1*v1.viscous_force_coefficients+v2.viscous_force_coefficients;
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void COUPLED_SYSTEM_VECTOR<TV>::
Print() const
{
    // Flat print
    for(int i=0;i<pressure.m;i++)
        LOG::cout<<pressure(i)<<" ";
    for(COUPLING_CONSTRAINT_ID i(0);i<lambda.Size();i++)
        LOG::cout<<lambda(i)<<" ";
    for(FORCE_AGGREGATE_ID i(0);i<force_coefficients.Size();i++)
        LOG::cout<<force_coefficients(i)<<" ";
    for(VISCOUS_FORCE_ID i(0);i<viscous_force_coefficients.Size();i++)
        LOG::cout<<viscous_force_coefficients(i)<<" ";
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int COUPLED_SYSTEM_VECTOR<TV>::
Raw_Size() const
{
    return pressure.m+Value(force_coefficients.m)+Value(lambda.m)+Value(viscous_force_coefficients.m);
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& COUPLED_SYSTEM_VECTOR<TV>::
Raw_Get(int i)
{
    if(i<pressure.m) return pressure(i);
    i-=pressure.m;
    int l=Value(lambda.m);
    if(i<l) return lambda(COUPLING_CONSTRAINT_ID(i));
    int f=Value(force_coefficients.m);
    if(i<l+f) return force_coefficients(FORCE_AGGREGATE_ID(i-l));
    return viscous_force_coefficients(VISCOUS_FORCE_ID(i-l-f));
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* COUPLED_SYSTEM_VECTOR<TV>::
Clone_Default() const
{
    COUPLED_SYSTEM_VECTOR<TV>* v=new COUPLED_SYSTEM_VECTOR<TV>;
    v->pressure.Resize(pressure.m);
    v->lambda.Resize(lambda.m);
    v->force_coefficients.Resize(force_coefficients.m);
    v->viscous_force_coefficients.Resize(viscous_force_coefficients.m);
    return v;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void COUPLED_SYSTEM_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& v)
{
    const COUPLED_SYSTEM_VECTOR<TV>& cs=debug_cast<const COUPLED_SYSTEM_VECTOR<TV>&>(v);
    pressure.Resize(cs.pressure.m);
    lambda.Resize(cs.lambda.m);
    force_coefficients.Resize(cs.force_coefficients.m);
    viscous_force_coefficients.Resize(cs.viscous_force_coefficients.m);
}
//#####################################################################
namespace PhysBAM{
template class COUPLED_SYSTEM_VECTOR<VECTOR<float,1> >;
template class COUPLED_SYSTEM_VECTOR<VECTOR<float,2> >;
template class COUPLED_SYSTEM_VECTOR<VECTOR<float,3> >;
template class COUPLED_SYSTEM_VECTOR<VECTOR<double,1> >;
template class COUPLED_SYSTEM_VECTOR<VECTOR<double,2> >;
template class COUPLED_SYSTEM_VECTOR<VECTOR<double,3> >;
}
