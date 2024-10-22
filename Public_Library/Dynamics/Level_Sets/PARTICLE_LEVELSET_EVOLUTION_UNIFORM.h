//#####################################################################
// Copyright 2003-2007, Ronald Fedkiw, Eran Guendelman, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_EVOLUTION_UNIFORM
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_EVOLUTION_UNIFORM__
#define __PARTICLE_LEVELSET_EVOLUTION_UNIFORM__

#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class TV> class RUNGEKUTTA;
template<class TV> class LEVELSET_ADVECTION;

template<class TV>
class PARTICLE_LEVELSET_EVOLUTION_UNIFORM:public PARTICLE_LEVELSET_EVOLUTION<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
    typedef ARRAY<RUNGEKUTTA<ARRAY_VIEW<TV> >*,TV_INT> T_ARRAYS_RUNGEKUTTA;
public:
    typedef PARTICLE_LEVELSET_EVOLUTION<T> BASE;
    using BASE::track_mass;using BASE::initial_mass;using BASE::runge_kutta_order_levelset;using BASE::runge_kutta_order_particles;using BASE::use_particle_levelset;
    using BASE::use_frozen_velocity;using BASE::cfl_number;using BASE::reseeding_frequency;using BASE::time;using BASE::use_fmm;using BASE::use_reinitialization;

    GRID<TV> grid;
    ARRAY<T,TV_INT> phi;
    ARRAY<T,FACE_INDEX<TV::m> > V;

private:
    PARTICLE_LEVELSET_UNIFORM<TV>* particle_levelset;
    LEVELSET_ADVECTION<TV>* levelset_advection;
public:

    PARTICLE_LEVELSET_EVOLUTION_UNIFORM(const GRID<TV>& grid_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,const int number_of_ghost_cells_input,bool multiphase);
    virtual ~PARTICLE_LEVELSET_EVOLUTION_UNIFORM();

    virtual PARTICLE_LEVELSET_UNIFORM<TV>& Particle_Levelset(const int i)
    {assert(i==0);return *particle_levelset;}

    virtual LEVELSET<TV>& Levelset(const int i)
    {assert(i==0);return particle_levelset->levelset;}
    
    virtual LEVELSET_ADVECTION<TV>& Levelset_Advection(const int i)
    {assert(i==0);return *levelset_advection;}

    virtual void Use_Semi_Lagrangian_Advection()
    {levelset_advection->Use_Local_Semi_Lagrangian_Advection();}
   
    virtual void Use_Hamilton_Jacobi_Weno_Advection()
    {levelset_advection->Use_Local_WENO_For_Advection();}

    virtual void Use_Hamilton_Jacobi_Eno_Advection(const int order)
    {assert(order>=1 && order<=3);levelset_advection->Use_Local_ENO_For_Advection(order);}

    virtual  void Track_Mass(const bool track_mass_input=true)
    {track_mass=track_mass_input;if(track_mass){initial_mass=levelset_advection->Approximate_Negative_Material();LOG::cout << "negative volume = " << initial_mass << std::endl;}}

    virtual void Initialize_Domain(const GRID<TV>& grid_input)
    {assert(grid_input.Is_MAC_Grid());grid=grid_input;phi.Resize(grid.Domain_Indices(particle_levelset->number_of_ghost_cells));V.Resize(grid);
    particle_levelset->Initialize_Particle_Levelset_Grid_Values();
    if(levelset_advection->semi_lagrangian_collidable) particle_levelset->levelset.Initialize_Valid_Masks(grid);}

    virtual void Make_Signed_Distance()
    {if(use_fmm){
            particle_levelset->levelset.Fast_Marching_Method(time,particle_levelset->levelset.half_band_width+grid.dX.Max()*(1+min(3,levelset_advection->local_advection_spatial_order)));
            particle_levelset->levelset.boundary->Apply_Boundary_Condition(grid,phi,time);}
    else if(use_reinitialization) levelset_advection->Reinitialize();}

    virtual void Set_Number_Particles_Per_Cell(const int number_particles_per_cell)
    {particle_levelset->Set_Number_Particles_Per_Cell(number_particles_per_cell);}

    virtual void Set_Levelset_Callbacks(LEVELSET_CALLBACKS<TV>& levelset_callbacks)
    {particle_levelset->levelset.Set_Levelset_Callbacks(levelset_callbacks);}

    virtual void Bias_Towards_Negative_Particles(const bool bias_towards_negative_particles)
    {particle_levelset->Bias_Towards_Negative_Particles(bias_towards_negative_particles);}

    virtual void Set_Seed(const int seed)
    {particle_levelset->random.Set_Seed(seed);}

    virtual void Seed_Particles(const T time)
    {particle_levelset->Seed_Particles(time);}

    virtual void Delete_Particles_Outside_Grid()
    {particle_levelset->Delete_Particles_Outside_Grid();}

    virtual void Fill_Levelset_Ghost_Cells(const T time)
    {particle_levelset->levelset.boundary->Fill_Ghost_Cells(grid,phi,particle_levelset->levelset.phi,0,time,particle_levelset->number_of_ghost_cells);}

    void Set_CFL_Number(const T cfl_number_input) override
    {PARTICLE_LEVELSET_EVOLUTION<T>::Set_CFL_Number(cfl_number_input);
    particle_levelset->cfl_number=cfl_number_input;}

//#####################################################################
    virtual void Advance_To_Time(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities,const T stopping_time,const bool verbose=true);
    virtual T Time_Step(const T stopping_time,bool& limited_by_stopping_time);
    virtual T CFL(const bool need_to_get_velocity=true,const bool analytic_test=false);
    virtual void Advance_One_Time_Step(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities,const T dt);
    virtual void Advance_Levelset(const T dt);
    virtual void Advance_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const bool analytic_test=false);
    virtual T Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time);
    virtual T Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time);
    virtual void Modify_Levelset_And_Particles(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities);
    virtual void Reseed_Particles(const T time,const int time_step=0,ARRAY<bool,TV_INT>* cell_centered_mask=0,const bool verbose=true);
//#####################################################################
};
}
#endif
