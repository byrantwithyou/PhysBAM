//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Forces_And_Torques/RIGID_GRAVITY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_GRAVITY<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    if(!gravity.Magnitude_Squared()) return;
    for(int p:force_rigid_body_particles){rigid_F(p).linear+=rigid_body_collection.rigid_body_particles.mass(p)*gravity;}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_GRAVITY<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(int p:force_rigid_body_particles){
        potential_energy-=rigid_body_collection.rigid_body_particles.mass(p)*TV::Dot_Product(rigid_body_collection.rigid_body_particles.frame(p).t,gravity);}
    return potential_energy;
}
//#####################################################################
namespace PhysBAM{
template class RIGID_GRAVITY<VECTOR<float,1> >;
template class RIGID_GRAVITY<VECTOR<float,2> >;
template class RIGID_GRAVITY<VECTOR<float,3> >;
template class RIGID_GRAVITY<VECTOR<double,1> >;
template class RIGID_GRAVITY<VECTOR<double,2> >;
template class RIGID_GRAVITY<VECTOR<double,3> >;
}
