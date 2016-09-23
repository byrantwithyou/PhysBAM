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
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();F(p)+=particles.mass(p)*gravity;}
    for(ELEMENT_ITERATOR iterator(force_rigid_body_particles);iterator.Valid();iterator.Next()){
        int p=iterator.Data();rigid_F(p).linear+=rigid_body_collection.rigid_body_particles.mass(p)*gravity;}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR GRAVITY<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        potential_energy-=particles.mass(p)*TV::Dot_Product(particles.X(p),gravity);}
    for(ELEMENT_ITERATOR iterator(force_rigid_body_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        potential_energy-=rigid_body_collection.rigid_body_particles.mass(p)*TV::Dot_Product(rigid_body_collection.rigid_body_particles.frame(p).t,gravity);}
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
