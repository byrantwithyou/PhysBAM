//#####################################################################
// Copyright 2008.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Dynamics/Coupled_Evolution/PRESSURE_VELOCITY_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PRESSURE_VELOCITY_VECTOR<TV>::
PRESSURE_VELOCITY_VECTOR(GENERALIZED_VELOCITY<TV>& solid_velocity_input,ARRAY<ARRAY<T> >& pressure_input)
    :solid_velocity(solid_velocity_input),pressure(pressure_input),deep_copy(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PRESSURE_VELOCITY_VECTOR<TV>::
~PRESSURE_VELOCITY_VECTOR()
{
    if(deep_copy){
        delete &solid_velocity;
        delete &pressure;}
}
//#####################################################################
// Function operator=
//#####################################################################
template<class TV> PRESSURE_VELOCITY_VECTOR<TV>& PRESSURE_VELOCITY_VECTOR<TV>::
operator=(const PRESSURE_VELOCITY_VECTOR& v)
{
    solid_velocity=v.solid_velocity;
    pressure=v.pressure;
    return *this;
}
//#####################################################################
// Function operator+=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& PRESSURE_VELOCITY_VECTOR<TV>::
operator+=(const BASE& bv)
{
    const PRESSURE_VELOCITY_VECTOR& v=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv);
    solid_velocity+=v.solid_velocity;
    pressure+=v.pressure;
    return *this;
}
//#####################################################################
// Function operator-=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& PRESSURE_VELOCITY_VECTOR<TV>::
operator-=(const BASE& bv)
{
    const PRESSURE_VELOCITY_VECTOR& v=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv);
    solid_velocity-=v.solid_velocity;
    pressure-=v.pressure;
    return *this;
}
//#####################################################################
// Function operator*=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& PRESSURE_VELOCITY_VECTOR<TV>::
operator*=(const T a)
{
    solid_velocity*=a;
    pressure*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void PRESSURE_VELOCITY_VECTOR<TV>::
Copy(const T c,const BASE& bv)
{
    const PRESSURE_VELOCITY_VECTOR& v=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv);
    assert(v.pressure.m==pressure.m);
    solid_velocity.Copy(c,v.solid_velocity);
    for(int i=0;i<v.pressure.m;i++) pressure(i).Copy(c,v.pressure(i));
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void PRESSURE_VELOCITY_VECTOR<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const PRESSURE_VELOCITY_VECTOR& v1=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv1),&v2=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv2);
    assert(v1.pressure.m==v2.pressure.m && pressure.m==v1.pressure.m);
    solid_velocity.Copy(c1,v1.solid_velocity,v2.solid_velocity);
    for(int i=0;i<v1.pressure.m;i++) pressure(i).Copy(c1,v1.pressure(i),v2.pressure(i));
}
//#####################################################################
// Function Dot
//#####################################################################
template<class TV> typename TV::SCALAR PRESSURE_VELOCITY_VECTOR<TV>::
Dot(const BASE& bv) const
{
    const PRESSURE_VELOCITY_VECTOR& v=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv);
    assert(pressure.m==v.pressure.m);
    T x=solid_velocity.Dot(v.solid_velocity);
    for(int i=0;i<v.pressure.m;i++) x+=pressure(i).Dot(v.pressure(i));
    return x;
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int PRESSURE_VELOCITY_VECTOR<TV>::
Raw_Size() const
{
    int n=solid_velocity.Raw_Size();
    for(int i=0;i<pressure.m;i++) n+=pressure(i).m;
    return n;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& PRESSURE_VELOCITY_VECTOR<TV>::
Raw_Get(int i)
{
    int n=solid_velocity.Raw_Size();
    if(i<n) return solid_velocity.Raw_Get(i);i-=n;
    for(int j=0;j<pressure.m;j++){
        if(i<pressure(j).m) return pressure(j)(i);
        i-=pressure(j).m;}
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>* PRESSURE_VELOCITY_VECTOR<TV>::
Clone_Default() const
{
    GENERALIZED_VELOCITY<TV>& gv=debug_cast<GENERALIZED_VELOCITY<TV>&>(*solid_velocity.Clone_Default());
    PRESSURE_VELOCITY_VECTOR<TV>* v=new PRESSURE_VELOCITY_VECTOR<TV>(gv,*new ARRAY<ARRAY<T> >(pressure.m));
    for(int i=0;i<v->pressure.m;i++) v->pressure(i).Resize(pressure(i).Size());
    v->deep_copy=true;
    return v;
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void PRESSURE_VELOCITY_VECTOR<TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& v)
{
    if(!deep_copy) return;
    const PRESSURE_VELOCITY_VECTOR<TV>& pv=debug_cast<const PRESSURE_VELOCITY_VECTOR<TV>&>(v);
    solid_velocity.Resize(pv.solid_velocity);
    pressure.Resize(pv.pressure.m);
    for(int i=0;i<pressure.m;i++) pressure(i).Resize(pv.pressure(i).Size());
}
namespace PhysBAM{
template class PRESSURE_VELOCITY_VECTOR<VECTOR<float,1> >;
template class PRESSURE_VELOCITY_VECTOR<VECTOR<float,2> >;
template class PRESSURE_VELOCITY_VECTOR<VECTOR<float,3> >;
template class PRESSURE_VELOCITY_VECTOR<VECTOR<double,1> >;
template class PRESSURE_VELOCITY_VECTOR<VECTOR<double,2> >;
template class PRESSURE_VELOCITY_VECTOR<VECTOR<double,3> >;
}
