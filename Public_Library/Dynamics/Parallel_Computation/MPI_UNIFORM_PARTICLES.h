//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_UNIFORM_PARTICLES
//#####################################################################
#ifndef __MPI_UNIFORM_PARTICLES__
#define __MPI_UNIFORM_PARTICLES__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>

template<class TV> class PARTICLE_LEVELSET;

namespace PhysBAM{

template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES>
void Exchange_Boundary_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset);

template<class TV,class T_PARTICLES>
void Exchange_Boundary_Particles_Flat(const MPI_UNIFORM_GRID<TV>& mpi_grid,T_PARTICLES& particles,const typename TV::SCALAR ghost_distance=0);

template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES>
void Exchange_Overlapping_Block_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset);

template<class TV,class T_PARTICLES,class T_ARRAYS_PARTICLES>
void Exchange_Ghost_Particles(const MPI_UNIFORM_GRID<TV>& mpi_grid,const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const int bandwidth,PARTICLE_LEVELSET<TV>& particle_levelset);

}
#endif
