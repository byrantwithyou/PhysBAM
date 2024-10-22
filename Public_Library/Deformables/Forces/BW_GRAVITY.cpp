//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Deformables/Forces/BW_GRAVITY.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void BW_GRAVITY<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(!gravity.Magnitude_Squared()) return;
    for(int p:force_particles) F(p)+=particles.mass(p)*gravity;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BW_GRAVITY<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(int p:force_particles){
        potential_energy-=particles.mass(p)*TV::Dot_Product(particles.X(p),gravity);}
    return potential_energy;
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
//TODO: This needs to be fixed if particles are deleted
template<class TV> void BW_GRAVITY<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    Update_Force_Particles(force_particles,IDENTITY_ARRAY<>(particles.Size()),particle_is_simulated,false);
}
//#####################################################################
namespace PhysBAM{
template class BW_GRAVITY<VECTOR<float,1> >;
template class BW_GRAVITY<VECTOR<float,2> >;
template class BW_GRAVITY<VECTOR<float,3> >;
template class BW_GRAVITY<VECTOR<double,1> >;
template class BW_GRAVITY<VECTOR<double,2> >;
template class BW_GRAVITY<VECTOR<double,3> >;
}
