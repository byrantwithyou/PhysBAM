//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void GRAVITY<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    if(!gravity.Magnitude_Squared()) return;
    for(int p:force_particles) F(p)+=particles.mass(p)*gravity;
    auto& rm=rigid_body_collection.rigid_body_particles.mass;
    for(int p:force_rigid_body_particles)
        rigid_F(p).linear+=rm(p)*gravity;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR GRAVITY<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(int p:force_particles)
        potential_energy-=particles.mass(p)*particles.X(p).Dot(gravity);
    auto& rbp=rigid_body_collection.rigid_body_particles;
    for(int p:force_rigid_body_particles)
        potential_energy-=rbp.mass(p)*rbp.frame(p).t.Dot(gravity);
    return potential_energy;
}
//#####################################################################
namespace PhysBAM{
template class GRAVITY<VECTOR<float,1> >;
template class GRAVITY<VECTOR<float,2> >;
template class GRAVITY<VECTOR<float,3> >;
template class GRAVITY<VECTOR<double,1> >;
template class GRAVITY<VECTOR<double,2> >;
template class GRAVITY<VECTOR<double,3> >;
}
