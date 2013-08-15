//#####################################################################
// Copyright 2007-2008, Michael Lentine, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/IDENTITY_ARRAY.h>
#include <Deformables/Forces/POINTWISE_DEFORMABLE_FORCE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Update_Mpi
//#####################################################################
//TODO: This needs to be fixed if particles are deleted
template<class TV> void POINTWISE_DEFORMABLE_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    this->mpi_solids=mpi_solids;
    if(influence_all_particles) force_particles.Update(Get_Particle_List(IDENTITY_ARRAY<>(particles.Size())),particle_is_simulated);
    else if(influenced_particles) force_particles.Update(Get_Particle_List(*influenced_particles),particle_is_simulated);
    is_simulated=particle_is_simulated;
}
//#####################################################################
namespace PhysBAM{
template class POINTWISE_DEFORMABLE_FORCE<VECTOR<float,1> >;
template class POINTWISE_DEFORMABLE_FORCE<VECTOR<float,2> >;
template class POINTWISE_DEFORMABLE_FORCE<VECTOR<float,3> >;
template class POINTWISE_DEFORMABLE_FORCE<VECTOR<double,1> >;
template class POINTWISE_DEFORMABLE_FORCE<VECTOR<double,2> >;
template class POINTWISE_DEFORMABLE_FORCE<VECTOR<double,3> >;
}
