//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VORTEX_PARTICLE_EVOLUTION_3D__
#define __VORTEX_PARTICLE_EVOLUTION_3D__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Incompressible/Incompressible_Flows/SCATTERED_INTERPOLATION.h>
#include <Incompressible/Particles/VORTICITY_PARTICLES.h>
namespace PhysBAM{
template<class T> class MPI_UNIFORM_GRID;

template<class T>
class VORTEX_PARTICLE_EVOLUTION_3D
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef MPI_UNIFORM_GRID<TV> T_MPI_GRID;
public:
    VORTICITY_PARTICLES<VECTOR<T,3> > vorticity_particles;
    GRID<TV> grid;
    T_MPI_GRID* mpi_grid;
    RANDOM_NUMBERS<T> random;
    SCATTERED_INTERPOLATION<TV> scattered_interpolation;
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > grid_vorticity,grid_vorticity_particles;
    T grid_confinement_parameter,particle_confinement_parameter,force_scaling;
    bool renormalize_vorticity_after_stretching_tilting;
    bool remove_grid_vorticity_from_particle_vorticity;
    bool apply_individual_particle_forces;
private:
//    T radius,radius_squared,one_over_radius_squared,one_over_radius_cubed;
public:

    VORTEX_PARTICLE_EVOLUTION_3D()
        :mpi_grid(0),grid_confinement_parameter(0),particle_confinement_parameter(0),force_scaling((T).03),renormalize_vorticity_after_stretching_tilting(false),
        remove_grid_vorticity_from_particle_vorticity(false),apply_individual_particle_forces(true)
    {
        scattered_interpolation.Use_Tent_Weights();
//        Set_Radius();
    }

    VORTEX_PARTICLE_EVOLUTION_3D(const VORTEX_PARTICLE_EVOLUTION_3D&) = delete;
    void operator=(const VORTEX_PARTICLE_EVOLUTION_3D&) = delete;

    void Initialize(const GRID<TV>& grid_input)
    {grid=grid_input.Get_MAC_Grid();grid_vorticity.Resize(grid.Domain_Indices(2),no_init);grid_vorticity_particles.Resize(grid.Domain_Indices(2),no_init);
    /*scattered_interpolation.Set_Radius_Of_Influence(radius);*/}

//#####################################################################
private:
    T Gaussian_Kernel(T distance_squared);
public:
//    void Set_Radius(const T radius_input=(T).01);
    void Compute_Body_Force(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& force,const T dt,const T time);
    void Compute_Body_Force(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time);
    void Euler_Step(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const T dt,const T time);
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const;
    void Read_Output_Files(const std::string& input_directory,const int frame);
//#####################################################################
};
}
#endif
