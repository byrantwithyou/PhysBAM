//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
using namespace PhysBAM;
template<class TV> ETHER_DRAG<TV>::
ETHER_DRAG(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_particles_input,
    ARRAY<int>* influenced_rigid_body_particles_input,T dynamic_ether_viscosity,T angular_viscosity)
    :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,influenced_particles_input,influenced_rigid_body_particles_input),use_constant_wind(dynamic_ether_viscosity!=0),
    constant_wind_viscosity(dynamic_ether_viscosity),constant_wind_angular_viscosity(angular_viscosity),use_spatially_varying_wind(false),spatially_varying_wind_viscosity(0)
{
}
template<class TV> ETHER_DRAG<TV>::
ETHER_DRAG(DEFORMABLE_PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_particles_input,
    const bool influence_all_rigid_body_particles_input,T dynamic_ether_viscosity,T angular_viscosity)
    :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,influence_all_particles_input,influence_all_rigid_body_particles_input),use_constant_wind(dynamic_ether_viscosity!=0),
    constant_wind_viscosity(dynamic_ether_viscosity),constant_wind_angular_viscosity(angular_viscosity),use_spatially_varying_wind(false),spatially_varying_wind_viscosity(0)
{
}
template<class TV> ETHER_DRAG<TV>::
~ETHER_DRAG()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void ETHER_DRAG<TV>::
Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    for(int k:force_particles){
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(particles.X(k))) F.V.array(k)+=spatially_varying_wind_viscosity*particles.mass(k)*Spatially_Varying_Wind_Velocity(particles.X(k));
            else if(use_constant_wind) F.V.array(k)+=constant_wind_viscosity*particles.mass(k)*constant_wind;}
        else if(use_constant_wind) F.V.array(k)+=constant_wind_viscosity*particles.mass(k)*constant_wind;}
    for(int k:force_rigid_body_particles){
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(rigid_body_collection.rigid_body_particles.frame(k).t))
                F.rigid_V.array(k).linear+=spatially_varying_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*
                    Spatially_Varying_Wind_Velocity(rigid_body_collection.rigid_body_particles.frame(k).t);
            else if(use_constant_wind) F.rigid_V.array(k).linear+=constant_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*constant_wind;}
        else if(use_constant_wind) F.rigid_V.array(k).linear+=constant_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*constant_wind;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void ETHER_DRAG<TV>::
Add_Velocity_Dependent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    for(int k:force_particles){
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(particles.X(k))) F.V.array(k)-=spatially_varying_wind_viscosity*particles.mass(k)*V.V.array(k);
            else if(use_constant_wind) F.V.array(k)-=constant_wind_viscosity*particles.mass(k)*V.V.array(k);}
        else if(use_constant_wind) F.V.array(k)-=constant_wind_viscosity*particles.mass(k)*V.V.array(k);}
    for(int k:force_rigid_body_particles){
        if(use_spatially_varying_wind){
            if(spatially_varying_wind_domain.Lazy_Inside(rigid_body_collection.rigid_body_particles.frame(k).t))
                F.rigid_V.array(k).linear-=spatially_varying_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*V.rigid_V.array(k).linear;
            else if(use_constant_wind) F.rigid_V.array(k).linear-=constant_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*V.rigid_V.array(k).linear;}
        else if(use_constant_wind) F.rigid_V.array(k).linear-=constant_wind_viscosity*rigid_body_collection.rigid_body_particles.mass(k)*V.rigid_V.array(k).linear;
        if(constant_wind_angular_viscosity) F.rigid_V.array(k).angular-=constant_wind_angular_viscosity*V.rigid_V.array(k).angular;}
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void ETHER_DRAG<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void ETHER_DRAG<TV>::
Add_Implicit_Velocity_Independent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time,bool transpose) const
{
}
//#####################################################################
namespace PhysBAM{
template class ETHER_DRAG<VECTOR<float,2> >;
template class ETHER_DRAG<VECTOR<float,3> >;
template class ETHER_DRAG<VECTOR<double,2> >;
template class ETHER_DRAG<VECTOR<double,3> >;
}
