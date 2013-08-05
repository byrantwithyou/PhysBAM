//#####################################################################
// Copyright 2005-2007, Ron Fedkiw, Frank Losasso, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_MULTIPHASE_UNIFORM
//#####################################################################
#ifndef __INCOMPRESSIBLE_MULTIPHASE_UNIFORM__
#define __INCOMPRESSIBLE_MULTIPHASE_UNIFORM__

#include <Incompressible/Incompressible_Flows/FLUID_STRAIN_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
namespace PhysBAM{

template<class TV>
class INCOMPRESSIBLE_MULTIPHASE_UNIFORM:public INCOMPRESSIBLE_UNIFORM<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef INTERPOLATION_UNIFORM<TV,T> T_INTERPOLATION_SCALAR;
    typedef ARRAY<typename TV::SPIN,TV_INT> T_ARRAYS_SPIN;
public:
    typedef INCOMPRESSIBLE_UNIFORM<TV> BASE;
    using BASE::use_force;using BASE::viscosity;using BASE::use_variable_viscosity;using BASE::use_variable_vorticity_confinement;using BASE::dt_old;using BASE::gravity;
    using BASE::downward_direction;using BASE::vorticity_confinements;using BASE::nonzero_viscosity;using BASE::nonzero_surface_tension;using BASE::mpi_grid;
    using BASE::use_explicit_part_of_implicit_viscosity;using BASE::vorticity_confinement;using BASE::max_time_step;using BASE::advection;
    using BASE::Set_Custom_Advection;using BASE::GFM;using BASE::number_of_interface_cells;using BASE::viscosities;using BASE::surface_tensions;
    using BASE::projection;using BASE::grid;using BASE::boundary;using BASE::force;using BASE::variable_vorticity_confinement;using BASE::strain;using BASE::variable_viscosity;
    using BASE::maximum_implicit_viscosity_iterations;using BASE::Extrapolate_Velocity_Across_Interface;

    ARRAY<T,FACE_INDEX<TV::m> > viscous_force;
    LEVELSET<TV>* levelset_for_dirichlet_regions;
    ARRAY<FLUID_STRAIN_UNIFORM<TV>*> strains;

    INCOMPRESSIBLE_MULTIPHASE_UNIFORM(const GRID<TV>& grid_input,PROJECTION_DYNAMICS_UNIFORM<TV>& projection_input);
    virtual ~INCOMPRESSIBLE_MULTIPHASE_UNIFORM();

    void Use_Strain(const ARRAY<bool>& use_multiphase_strain)
    {for(int i=0;i<strains.m;i++)delete strains(i);
    strains.Resize(use_multiphase_strain.m);
    for(int i=0;i<use_multiphase_strain.m;i++)if(use_multiphase_strain(i))strains(i)=new FLUID_STRAIN_UNIFORM<TV>(grid);}

    // overrides version from BASE
    void Advance_One_Time_Step_Forces(const T dt,const T time,const bool implicit_viscosity=false,const ARRAY<T,TV_INT>* phi_ghost=0)
    {PHYSBAM_NOT_IMPLEMENTED();/*PHYSBAM_ASSERT(!phi_ghost);Advance_One_Time_Step_Forces(dt,time,implicit_viscosity,0,0);*/}

    // overrides version from BASE
    void Advance_One_Time_Step_Convection(const T dt,const T time,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_to_advect)
    {PHYSBAM_NOT_IMPLEMENTED();/*Advance_One_Time_Step_Convection(dt,time,face_velocities_to_advect,0);*/}
    
//#####################################################################
    void Advance_One_Time_Step_Forces(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,const bool implicit_viscosity,const ARRAY<ARRAY<T,TV_INT>>* phi_ghost,
        const ARRAY<bool>* pseudo_dirichlet_regions,const int number_of_ghost_cells);
    void Advance_One_Time_Step_Convection(const T dt,const T time,ARRAY<T,FACE_INDEX<TV::m> >& advecting_face_velocities,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_to_advect,const ARRAY<bool>* pseudo_dirichlet_regions,const int number_of_ghost_cells);
    void Advance_One_Time_Step_Implicit_Part(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,const bool implicit_viscosity=false);
    void Calculate_Pressure_Jump(const T dt,const T time);
    T CFL(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool inviscid=false,const bool viscous_only=false) const;
    void Set_Dirichlet_Boundary_Conditions(ARRAY<ARRAY<T,TV_INT>>& phis,const ARRAY<bool>& dirichlet_regions,const ARRAY<T>* pressures=0);
    void Add_Surface_Tension(LEVELSET<TV>& levelset,const T time);
    void Compute_Vorticity_Confinement_Force(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<TV,TV_INT>& F) PHYSBAM_OVERRIDE;
protected:
    void Discretize_Explicit_Viscous_Terms(const T dt){PHYSBAM_NOT_IMPLEMENTED();}
    void Implicit_Viscous_Update(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
//#####################################################################
};
}
#endif
