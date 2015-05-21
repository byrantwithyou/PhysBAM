//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Forces/ELASTIC_ETHER_DRAG.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void ELASTIC_ETHER_DRAG<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){
        int k=iterator.Data();
        F(k)-=coefficient*particles.mass(k)*particles.V(k);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void ELASTIC_ETHER_DRAG<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void ELASTIC_ETHER_DRAG<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void ELASTIC_ETHER_DRAG<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const
{
    T c=coefficient*dt_dv_over_dx/dt*scale;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){
        int k=iterator.Data();
        F(k)-=c*particles.mass(k)*V(k);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void ELASTIC_ETHER_DRAG<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR ELASTIC_ETHER_DRAG<TV>::
Potential_Energy(const T time) const
{
    T pe=0;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){
        int k=iterator.Data();
        pe+=particles.mass(k)*particles.V(k).Magnitude_Squared();}
    return pe*coefficient*dt/(2*dt_dv_over_dx);
}
//#####################################################################
namespace PhysBAM{
template class ELASTIC_ETHER_DRAG<VECTOR<float,2> >;
template class ELASTIC_ETHER_DRAG<VECTOR<float,3> >;
template class ELASTIC_ETHER_DRAG<VECTOR<double,2> >;
template class ELASTIC_ETHER_DRAG<VECTOR<double,3> >;
}
