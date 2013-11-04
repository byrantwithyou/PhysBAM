//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Forces_And_Torques/RIGID_ETHER_DRAG.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> RIGID_ETHER_DRAG<TV>::
RIGID_ETHER_DRAG(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_rigid_body_particles_input,T dynamic_ether_viscosity,T angular_viscosity)
    :RIGID_POINTWISE_FORCE<TV>(rigid_body_collection_input,influenced_rigid_body_particles_input),use_constant_wind(dynamic_ether_viscosity!=0),
    constant_wind_viscosity(dynamic_ether_viscosity),constant_wind_angular_viscosity(angular_viscosity),use_spatially_varying_wind(false),spatially_varying_wind_viscosity(0)
{
}
template<class TV> RIGID_ETHER_DRAG<TV>::
RIGID_ETHER_DRAG(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_rigid_body_particles_input,T dynamic_ether_viscosity,T angular_viscosity)
    :RIGID_POINTWISE_FORCE<TV>(rigid_body_collection_input,influence_all_rigid_body_particles_input),use_constant_wind(dynamic_ether_viscosity!=0),
    constant_wind_viscosity(dynamic_ether_viscosity),constant_wind_angular_viscosity(angular_viscosity),use_spatially_varying_wind(false),spatially_varying_wind_viscosity(0)
{
}
template<class TV> RIGID_ETHER_DRAG<TV>::
~RIGID_ETHER_DRAG()
{
}
template<class TV> void RIGID_ETHER_DRAG<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_rigid_body_particles);iterator.Valid();iterator.Next()){int k=iterator.Data();
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(rigid_body_collection.rigid_body_particles.frame(k).t))
                rigid_F(k).linear+=spatially_varying_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*
                    Spatially_Varying_Wind_Velocity(rigid_body_collection.rigid_body_particles.frame(k).t);
            else if(use_constant_wind) rigid_F(k).linear+=constant_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*constant_wind;}
        else if(use_constant_wind) rigid_F(k).linear+=constant_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*constant_wind;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RIGID_ETHER_DRAG<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_rigid_body_particles);iterator.Valid();iterator.Next()){int k=iterator.Data();
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(rigid_body_collection.rigid_body_particles.frame(k).t))
                rigid_F(k).linear-=spatially_varying_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*rigid_V(k).linear;
            else if(use_constant_wind) rigid_F(k).linear-=constant_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*rigid_V(k).linear;}
        else if(use_constant_wind) rigid_F(k).linear-=constant_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*rigid_V(k).linear;
        if(constant_wind_angular_viscosity) rigid_F(k).angular-=constant_wind_angular_viscosity*rigid_V(k).angular;}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_ETHER_DRAG<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void RIGID_ETHER_DRAG<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_ETHER_DRAG<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
}
//#####################################################################
namespace PhysBAM{
template class RIGID_ETHER_DRAG<VECTOR<float,2> >;
template class RIGID_ETHER_DRAG<VECTOR<float,3> >;
template class RIGID_ETHER_DRAG<VECTOR<double,2> >;
template class RIGID_ETHER_DRAG<VECTOR<double,3> >;
}
