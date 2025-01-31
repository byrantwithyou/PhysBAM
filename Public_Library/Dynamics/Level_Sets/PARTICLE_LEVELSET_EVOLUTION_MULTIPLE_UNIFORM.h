//#####################################################################
// Copyright 2005-2006, Ron Fedkiw, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM__
#define __PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM__

#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_MULTIPLE_UNIFORM.h>
namespace PhysBAM{

template<class TV> class LEVELSET_ADVECTION_MULTIPLE;

template<class TV>
class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM:public PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
public:
    typedef PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV> BASE;
    using BASE::track_mass;using BASE::runge_kutta_order_levelset;using BASE::runge_kutta_order_particles;using BASE::use_particle_levelset;using BASE::use_frozen_velocity;
    using BASE::cfl_number;using BASE::reseeding_frequency;using BASE::time;using BASE::use_fmm;using BASE::use_reinitialization;using BASE::V;using BASE::grid;

    ARRAY<ARRAY<T,TV_INT>> phis;
    PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>& particle_levelset_multiple;
    ARRAY<T> initial_mass;

    LEVELSET_ADVECTION_MULTIPLE<TV>& levelset_advection_multiple;

    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM(const GRID<TV>& grid_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,const int number_of_ghost_cells_input);
    virtual ~PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM();

//#####################################################################
    virtual PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>& Particle_Levelset_Multiple();
    virtual LEVELSET_MULTIPLE<TV>& Levelset_Multiple();
    virtual PARTICLE_LEVELSET_UNIFORM<TV>& Particle_Levelset(const int i) override;
    virtual LEVELSET<TV>& Levelset(const int i) override;
    virtual LEVELSET_ADVECTION<TV>& Levelset_Advection(const int i) override;
    void Use_Semi_Lagrangian_Advection() override;
    void Use_Hamilton_Jacobi_Weno_Advection() override;
    void Use_Hamilton_Jacobi_Eno_Advection(const int order) override;
    void Track_Mass(const bool track_mass_input=true) override;
    void Initialize_Domain(const GRID<TV>& grid_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,const int number_of_regions,const bool use_only_negative_particles=true);  // don't call up to the base class here because we don't need those variables initialized OVERRIDE PROBLEM
    void Initialize_Domain(const GRID<TV>& grid_input) override;
    void Make_Signed_Distance() override;
    void Set_Number_Particles_Per_Cell(const int number_particles_per_cell) override;
    void Set_Levelset_Callbacks(LEVELSET_CALLBACKS<TV>& levelset_callbacks) override;
    void Bias_Towards_Negative_Particles(const bool bias_towards_negative_particles) override;
    void Set_Seed(const int seed) override;
    void Seed_Particles(const T time) override;
    void Delete_Particles_Outside_Grid() override;
    void Set_CFL_Number(const T cfl_number_input) override;
    void Advance_To_Time(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities,const T stopping_time,const bool verbose=true) override;
    T Time_Step(const T stopping_time,bool& limited_by_stopping_time) override;
    T CFL(const bool need_to_get_velocity=true,const bool analytic_test=false) override;
    void Advance_One_Time_Step(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities,const T dt) override;
    void Advance_Levelset(const T dt) override;
    void Advance_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const bool analytic_test=false) override;
    T Advance_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time) override;
    T Advance_Particles(ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time) override;
    void Modify_Levelset_And_Particles(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities) override;
    void Reseed_Particles(const T time,const int time_step=0,ARRAY<bool,TV_INT>* cell_centered_mask=0,const bool verbose=true) override;
    void Fill_Levelset_Ghost_Cells(const T time) override;
//#####################################################################
};
}
#endif
