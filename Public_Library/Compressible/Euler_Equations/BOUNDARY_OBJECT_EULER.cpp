//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Compressible/Euler_Equations/BOUNDARY_OBJECT_EULER.h>
#include <Compressible/Euler_Equations/EULER.h>

using namespace PhysBAM;
//#####################################################################
// Function Apply_Neumann_Boundary_Condition
//#####################################################################
template<class TV> void BOUNDARY_OBJECT_EULER<TV>::
Apply_Neumann_Boundary_Condition(TV_DIMENSION& u_1d,const T neumann_face_velocity,const int axis)
{
    T rho=u_1d(0);
    if(rho){
        TV velocity=EULER<TV>::Get_Velocity(u_1d);velocity(axis)=(T)2*neumann_face_velocity-velocity(axis);
        T internal_energy=EULER<TV>::e(u_1d);
        u_1d=EULER<TV>::Get_Euler_State_From_rho_velocity_And_internal_energy(rho,velocity,internal_energy);}
}
template<class TV> void BOUNDARY_OBJECT_EULER<TV>::
Apply_Neumann_Boundary_Condition(T_ARRAYS_ELEMENT& u_1d,const TV& normal,const T object_velocity_normal_component)
{
    T rho=u_1d(0);
    if(rho){
        TV velocity=EULER<TV>::Get_Velocity(u_1d);
        T velocity_normal_component=TV::Dot_Product(velocity,normal);
        TV velocity_tangential=velocity-normal*velocity_normal_component;
        velocity=velocity_tangential+normal*((T)2*object_velocity_normal_component-velocity_normal_component);
        T internal_energy=EULER<TV>::e(u_1d);
        u_1d=EULER<TV>::Get_Euler_State_From_rho_velocity_And_internal_energy(rho,velocity,internal_energy);}
}
//#####################################################################
namespace PhysBAM{
template class BOUNDARY_OBJECT_EULER<VECTOR<float,1> >;
template class BOUNDARY_OBJECT_EULER<VECTOR<float,2> >;
template class BOUNDARY_OBJECT_EULER<VECTOR<float,3> >;
template class BOUNDARY_OBJECT_EULER<VECTOR<double,1> >;
template class BOUNDARY_OBJECT_EULER<VECTOR<double,2> >;
template class BOUNDARY_OBJECT_EULER<VECTOR<double,3> >;
}
