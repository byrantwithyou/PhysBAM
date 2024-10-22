//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_DYNAMICS_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_DYNAMICS_UNIFORM__
#define __PROJECTION_DYNAMICS_UNIFORM__

#include <Core/Data_Structures/TRIPLE.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/PROJECTION_COLLIDABLE_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS.h>
namespace PhysBAM{

template<class TV> class DETONATION_SHOCK_DYNAMICS;
template<class TV> class FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM;

template<class TV>
class PROJECTION_DYNAMICS_UNIFORM:public PROJECTION_COLLIDABLE_UNIFORM<TV>,public PROJECTION_DYNAMICS<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;typedef FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<TV> T_FACE_LOOKUP_FIRE_MULTIPHASE;
public:
    typedef PROJECTION_COLLIDABLE_UNIFORM<TV> BASE;
    using BASE::use_non_zero_divergence;using BASE::p_grid;using BASE::poisson;using BASE::elliptic_solver;using BASE::laplace;using BASE::p;using BASE::collidable_solver;using BASE::use_divergence_multiplier;using BASE::divergence;
    using BASE::divergence_multiplier;using BASE::poisson_collidable;using BASE::laplace_collidable;
    using BASE::Use_Divergence_Multiplier;using BASE::Use_Non_Zero_Divergence;using BASE::Compute_Divergence;using BASE::density;
    using PROJECTION_DYNAMICS<T>::flame;using PROJECTION_DYNAMICS<T>::flame_speed_constants;using BASE::Enforce_Velocity_Compatibility;
    
    bool use_flame_speed_multiplier;
    ARRAY<T,FACE_INDEX<TV::m> > flame_speed_multiplier;
    DETONATION_SHOCK_DYNAMICS<TV>* dsd;

protected:
    bool use_divergence_multiplier_save_for_sph,use_non_zero_divergence_save_for_sph;
    ARRAY<T,TV_INT> *p_save_for_sph,*divergence_save_for_sph,*divergence_multiplier_save_for_sph;
    ARRAY<T,FACE_INDEX<TV::m> > *face_velocities_save_for_sph;
    LAPLACE_UNIFORM<TV>* elliptic_solver_save_for_sph;
    LAPLACE_COLLIDABLE_UNIFORM<TV>* laplace_save_for_sph;
    POISSON_COLLIDABLE_UNIFORM<TV>* poisson_save_for_sph;
    LAPLACE_COLLIDABLE<TV>* collidable_solver_save_for_sph;
public:

    PROJECTION_DYNAMICS_UNIFORM(const GRID<TV>& mac_grid,const bool flame_input=false,const bool multiphase=false,const bool use_variable_beta=false,const bool use_poisson=false);
    PROJECTION_DYNAMICS_UNIFORM(const GRID<TV>& mac_grid,LEVELSET<TV>& levelset_input);
    virtual ~PROJECTION_DYNAMICS_UNIFORM();

    T Face_Velocity_With_Ghost_Value_Multiphase(const T_ARRAYS_BASE& face_velocities_ghost,const int axis,const TV_INT& face_index,const int current_region) const
    {assert(flame);return face_velocities_ghost(face_index)-Face_Jump_Multiphase(axis,face_index,current_region);}

    T Face_Velocity_With_Ghost_Value_Multiphase(const T_ARRAYS_BASE& face_velocities_ghost,const int axis,const TV_INT& face_index,const int current_region,const int face_region) const
    {assert(flame);return face_velocities_ghost(face_index)-Face_Jump_Multiphase(axis,face_index,current_region,face_region);}

    T Face_Velocity_With_Ghost_Value_Multiphase(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const int axis,const TV_INT& face_index,const int current_region) const
    {return Face_Velocity_With_Ghost_Value_Multiphase(face_velocities_ghost.Component(axis),axis,face_index,current_region);}

    T Face_Velocity_With_Ghost_Value_Multiphase(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const int axis,const TV_INT& face_index,const int current_region,const int face_region) const
    {return Face_Velocity_With_Ghost_Value_Multiphase(face_velocities_ghost.Component(axis),axis,face_index,current_region,face_region);}

    T Face_Jump_Multiphase(const int axis,const TV_INT& face_index,const int current_region) const
    {int face_region=poisson_collidable->levelset_multiple->Inside_Region_Face(axis,face_index);return Face_Jump_Multiphase(axis,face_index,current_region,face_region);}

//#####################################################################
    virtual void Initialize_Grid(const GRID<TV>& mac_grid);
    void Initialize_Dsd(const LEVELSET_MULTIPLE<TV>& levelset_multiple,const ARRAY<bool>& is_fuel_region);
    void Initialize_Dsd(const LEVELSET<TV>& levelset,const ARRAY<bool>& fuel_region);
    virtual void Make_Divergence_Free(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    void Compute_Divergence(const T_FACE_LOOKUP_FIRE_MULTIPHASE& face_lookup,LAPLACE_UNIFORM<TV>* solver);
    void Apply_Pressure(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,bool scale_by_dt=false);
    void Set_Up_For_SPH(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool use_variable_density_solve=false,const bool use_one_way_coupling=false);
    void Restore_After_SPH(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool use_variable_density_solve=false,const bool use_one_way_coupling=false);
    void Update_Phi_And_Move_Velocity_Discontinuity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,LEVELSET_MULTIPLE<TV>& levelset_multiple,const T time,const bool update_phi_only=false);
    template<class FACE_LOOKUP> void Compute_Divergence(const FACE_LOOKUP &face_lookup,LAPLACE_UNIFORM<TV>* solver);
    T Flame_Speed_Face_Multiphase(const int axis,const TV_INT& face_index,const int fuel_region,const int product_region) const;
    void Use_Flame_Speed_Multiplier(const bool use_flame_speed_multiplier_input=true);
    T Face_Jump_Multiphase(const int axis,const TV_INT& face_index,const int current_region,const int face_region) const;
//#####################################################################
};
}
#endif
