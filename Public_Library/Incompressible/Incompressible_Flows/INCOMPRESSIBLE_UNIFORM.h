//#####################################################################
// Copyright 2002-2010, Mridul Aanjaneya, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Duc Nguyen, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_UNIFORM  
//#####################################################################
#ifndef __INCOMPRESSIBLE_UNIFORM__
#define __INCOMPRESSIBLE_UNIFORM__

#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_COLLIDABLE_POLICY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/FLUID_STRAIN_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
namespace PhysBAM{

template<class TV> class BOUNDARY_CONDITIONS_CALLBACKS;
template<class TV,class T> class EXTRAPOLATION_UNIFORM;
template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class TV>
class INCOMPRESSIBLE_UNIFORM:public INCOMPRESSIBLE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename ADVECTION_COLLIDABLE_POLICY<TV>::ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef EXTRAPOLATION_UNIFORM<TV,T> T_EXTRAPOLATION_SCALAR;typedef MPI_UNIFORM_GRID<TV> T_MPI_GRID;
    typedef typename TV::SPIN T_SPIN;typedef typename ARRAY<T,TV_INT>::template REBIND<T_SPIN>::TYPE T_ARRAYS_SPIN;
public:
    typedef INCOMPRESSIBLE<TV> BASE;
    using BASE::use_force;using BASE::surface_tension;using BASE::use_variable_surface_tension;using BASE::viscosity;using BASE::use_variable_viscosity;
    using BASE::use_variable_vorticity_confinement;using BASE::dt_old;using BASE::gravity;using BASE::nonzero_surface_tension;using BASE::nonzero_viscosity;
    using BASE::downward_direction;using BASE::Set_Custom_Advection;using BASE::valid_mask;
    using BASE::use_explicit_part_of_implicit_viscosity;using BASE::vorticity_confinement;using BASE::max_time_step;using BASE::advection;
    using BASE::maximum_implicit_viscosity_iterations;
    
    GRID<TV> grid;
    T_MPI_GRID* mpi_grid;
    BOUNDARY<TV,T>* boundary;
    PROJECTION_DYNAMICS_UNIFORM<TV>& projection;
    ARRAY<T,TV_INT> variable_surface_tension;
    ARRAY<T,TV_INT> variable_viscosity;
    ARRAY<T,FACE_INDEX<TV::m> > force;
    ARRAY<T,FACE_INDEX<TV::m> > potential_energy;
    ARRAY<T,TV_INT> variable_vorticity_confinement;
    FLUID_STRAIN_UNIFORM<TV>* strain;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_body_list;
    bool momentum_conserving_vorticity;    
    bool use_vorticity_weights;
    ARRAY<T,FACE_INDEX<TV::m> > vorticity_weights;
    T energy_clamp;
    int vc_projection_direction;
    T buoyancy_constant;
    THREAD_QUEUE* thread_queue;
protected:               
    BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>& boundary_default;
    ADVECTION_MACCORMACK_UNIFORM<TV,T,ADVECTION<TV,T> >* advection_maccormack;
public:

    INCOMPRESSIBLE_UNIFORM(const GRID<TV>& grid_input,PROJECTION_DYNAMICS_UNIFORM<TV>& projection_input,THREAD_QUEUE* thread_queue_input=0);
    virtual ~INCOMPRESSIBLE_UNIFORM();

    void Set_Custom_Boundary(BOUNDARY<TV,T>& boundary_input)
    {boundary=&boundary_input;}

    void Set_Collision_Body_List(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input)
    {collision_body_list=&collision_body_list_input;}

//#####################################################################
    void Initialize_Grids(const GRID<TV>& grid_input) PHYSBAM_OVERRIDE;
    void Set_Body_Force(const bool use_force_input=true);
    void Set_Variable_Surface_Tension(const bool use_variable_surface_tension_input=true);
    void Set_Variable_Viscosity(const bool use_variable_viscosity_input=true);
    void Use_Variable_Vorticity_Confinement(const bool use_variable_vorticity_confinement_input=true);
    void Use_Variable_Vorticity_Confinement(GRID<TV>& grid,const bool use_variable_vorticity_confinement_input=true);
    void Use_Strain();
    void Apply_Pressure_Kinetic_Energy(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_new,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_old,const T dt,const T time);
    void Add_Energy_With_Vorticity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const VECTOR<VECTOR<bool,2>,TV::dimension>& domain_boundary,const T dt,const T time,const int number_of_ghost_cells,LEVELSET<TV>* lsv=0,ARRAY<T,TV_INT>* density=0);
    void Advance_One_Time_Step_Convection(const T dt,const T time,const ARRAY<T,FACE_INDEX<TV::m> >& advecting_face_velocities,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_to_advect,const int number_of_ghost_cells);
    void Advance_One_Time_Step_Forces(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,const bool implicit_viscosity,const ARRAY<T,TV_INT>* phi_ghost,const int number_of_ghost_cells);
    void Add_Gravity_Threaded(RANGE<TV_INT>&domain,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,int axis);
    void Add_Body_Force_Threaded(RANGE<TV_INT>&domain,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,int axis);
    void Advance_One_Time_Step_Implicit_Part(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,const bool implicit_viscosity=false,BOUNDARY<TV,T>* projection_boundary=0,
        bool use_levelset_viscosity=false,BOUNDARY_CONDITIONS_CALLBACKS<TV>* bc_callbacks=0,bool print_viscosity_matrix=false);
    int Real_CFL(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool inviscid,const bool viscous_only,T input_dt) const;
    T CFL(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool inviscid=false,const bool viscous_only=false) const;
    void Apply_Vorticity_Confinement_Force(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<TV,TV_INT>& F);
    virtual void Compute_Vorticity_Confinement_Force(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<TV,TV_INT>& F);
    void Extrapolate_Velocity_Across_Interface(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_input,const ARRAY<T,TV_INT>& phi_ghost,const bool enforce_divergence_free=false,const T band_width=3,
        const T damping=0,const TV& air_speed=TV(),const ARRAY<VECTOR<bool,TV::m>,FACE_INDEX<TV::m> >* face_neighbors_visible=0,const ARRAY<bool,FACE_INDEX<TV::m> >* fixed_faces_input=0);
    void Set_Dirichlet_Boundary_Conditions(const ARRAY<T,TV_INT>* phi,const T pressure);
    void Set_Dirichlet_Boundary_Conditions(const ARRAY<T,TV_INT>* phi,const ARRAY<T,TV_INT>& pressure);
    void Add_Surface_Tension(LEVELSET<TV>& levelset,const T time);
    void Use_Maccormack_Advection(const ARRAY<bool,TV_INT>* node_mask, const ARRAY<bool,TV_INT>* cell_mask, const ARRAY<bool,FACE_INDEX<TV::m> >* face_mask);
    bool Consistent_Boundary_Conditions(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities) const;
    void Implicit_Viscous_Update(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
//#####################################################################
};
}
#endif

