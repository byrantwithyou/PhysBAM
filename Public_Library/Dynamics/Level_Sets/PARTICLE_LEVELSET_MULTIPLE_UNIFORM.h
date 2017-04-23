//#####################################################################
// Copyright 2005, Ron Fedkiw, Eran Guendelman, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_MULTIPLE_UNIFORM  
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_MULTIPLE_UNIFORM__
#define __PARTICLE_LEVELSET_MULTIPLE_UNIFORM__

#include <Incompressible/Level_Sets/LEVELSET_MULTIPLE.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class PARTICLE_LEVELSET_MULTIPLE_UNIFORM
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
public:
    T min_collision_distance_factor,max_collision_distance_factor,max_minus_min_collision_distance_factor_over_max_short;

    LEVELSET_MULTIPLE<TV> levelset_multiple;
    ARRAY<PARTICLE_LEVELSET_UNIFORM<TV>*> particle_levelsets;
    int number_of_ghost_cells;

    PARTICLE_LEVELSET_MULTIPLE_UNIFORM(GRID<TV>& grid_input,ARRAY<ARRAY<T,TV_INT>>& phis_input,const int number_of_ghost_cells_input);
    PARTICLE_LEVELSET_MULTIPLE_UNIFORM(const PARTICLE_LEVELSET_MULTIPLE_UNIFORM&) = delete;
    void operator=(const PARTICLE_LEVELSET_MULTIPLE_UNIFORM&) = delete;
    ~PARTICLE_LEVELSET_MULTIPLE_UNIFORM();

    void Initialize_Particle_Levelsets_And_Grid_Values(GRID<TV>& grid,ARRAY<ARRAY<T,TV_INT>>& phis,
        GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,
        const int number_of_regions,const bool only_use_negative_particles=true);
    void Store_Unique_Particle_Id();
    void Set_Band_Width(const T number_of_cells=6) ;
    void Use_Removed_Negative_Particles(const bool use_removed_negative_particles_input=true);
    void Use_Removed_Positive_Particles(const bool use_removed_positive_particles_input=true);
    void Set_Collision_Distance_Factors(const T min_collision_distance_factor_input=(T).1,const T max_collision_distance_factor_input=(T)1);
    T Particle_Collision_Distance(const unsigned short quantized_collision_distance);
    void Seed_Particles(const T time,const bool verbose=true);
    void Adjust_Particle_Radii();
    void Modify_Levelset_Using_Escaped_Particles(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities);
    void Euler_Step_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,const T dt,const T time=0,const bool use_second_order_for_nonremoved_particles=false,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    void Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    int Reseed_Particles(const T time,ARRAY<bool,TV_INT>* cell_centered_mask=0);
    void Delete_Particles_Outside_Grid();
    void Identify_And_Remove_Escaped_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,const T radius_fraction=1.5,const bool verbose=true);
    void Reincorporate_Removed_Particles(const T radius_fraction,const T mass_scaling,ARRAY<T,FACE_INDEX<TV::m> >* V,const bool conserve_momentum_for_removed_negative_particles=true);
//#####################################################################
};
}
#endif
