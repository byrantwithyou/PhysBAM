//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_UNIFORM
//#####################################################################
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_1D class.
// Input U as 3 by (1,m) for mass, momentum, and energy.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m).
// When psi=true, solve the equaitions.
// When psi=false, do NOT solve the equations.
//
//#####################################################################
#ifndef __EULER_UNIFORM__
#define __EULER_UNIFORM__

#include <Compressible/Euler_Equations/EULER.h>
#include <Compressible/Euler_Equations/EULER_CAVITATION_UNIFORM.h>
#include <Compressible/Euler_Equations/EULER_PROJECTION_UNIFORM.h>
namespace PhysBAM{
template<class T> class MPI_UNIFORM_GRID;

template<class TV>
class EULER_UNIFORM:public EULER<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;
    typedef ARRAY<TV_DIMENSION,FACE_INDEX<TV::m> > T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef MPI_UNIFORM_GRID<TV> T_MPI_GRID;
    typedef EULER<TV> BASE;
    typedef TV_DIMENSION T_ARRAYS_ELEMENT;
    typedef BOUNDARY<TV,TV_DIMENSION> T_BOUNDARY;
protected:
    using BASE::max_time_step;using BASE::cut_out_grid;using BASE::gravity;using BASE::downward_direction;
public:
    using BASE::boundary;using BASE::conservation;using BASE::eos;using BASE::Get_Velocity;using BASE::e;
    using BASE::Set_Max_Time_Step;using BASE::Set_Custom_Conservation;using BASE::Set_CFL_Number;using BASE::open_boundaries;
    using BASE::use_force;using BASE::Get_Density;

    GRID<TV> grid;
    T_MPI_GRID* mpi_grid;
    T_ARRAYS_DIMENSION_SCALAR U,U_save; // mass, momentum, and energy
    const T_ARRAYS_DIMENSION_SCALAR& U_ghost;
    ARRAY<bool,TV_INT>* psi_pointer; // defines cut out grid
    ARRAY<bool,TV_INT> psi;
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,TV::m> eigensystems,eigensystems_default,eigensystems_pressureonly;
    bool timesplit,use_sound_speed_for_cfl,perform_rungekutta_for_implicit_part,compute_pressure_fluxes,thinshell;
    bool use_sound_speed_based_dt_multiple_for_cfl; // if set, dt will be set to multiplication_factor_for_sound_speed_based_dt*dt_based_on_c whenever this number is less than dt_based_on_u
    T multiplication_factor_for_sound_speed_based_dt;
    EULER_PROJECTION_UNIFORM<TV> euler_projection;
    bool apply_cavitation_correction;
    EULER_CAVITATION_UNIFORM<TV> euler_cavitation_density;
    EULER_CAVITATION_UNIFORM<TV> euler_cavitation_internal_energy;
    T e_min,last_dt;
    TV_INT pressure_flux_dimension_indices;
    ARRAY<T,FACE_INDEX<TV::m> > force;
    TV_DIMENSION initial_total_conserved_quantity,accumulated_boundary_flux;
private:
    mutable T_ARRAYS_DIMENSION_SCALAR U_ghost_private; // Forces us to only touch U_ghost through the const reference and let Fill_Ghost_Cells be the only thing modifying the variable;
    T_FACE_ARRAYS_DIMENSION_SCALAR fluxes_pressure;
    ARRAY<T,TV_INT> added_internal_energy;
    mutable bool ghost_cells_valid;
    mutable int ghost_cells_valid_ring;
public:
    bool need_to_remove_added_internal_energy,need_to_remove_added_internal_energy_save;

    EULER_UNIFORM(const GRID<TV>& grid_input);
    virtual ~EULER_UNIFORM();

    T Get_Temperature(const TV_INT& cell)
    {
        T density=Get_Density(U,cell);
        T internal_energy=e(U,cell);
        return eos->T(density,internal_energy);
    }

//#####################################################################
    void Set_Up_Cut_Out_Grid(ARRAY<bool,TV_INT>& psi_input);
    void Set_Custom_Equation_Of_State(EOS<T>& eos_input);
    void Set_Custom_Boundary(T_BOUNDARY& boundary_input);
    void Set_Body_Force(const bool use_force_input=true);
    void Initialize_Domain(const GRID<TV>& grid_input);
    void Save_State(T_ARRAYS_DIMENSION_SCALAR& U_save,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_save,bool& need_to_remove_added_internal_energy_save);
    void Restore_State(T_ARRAYS_DIMENSION_SCALAR& U_save,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_save,bool& need_to_remove_added_internal_energy_save);
    void Get_Cell_Velocities(const T dt,const T time,const int ghost_cells,ARRAY<TV,TV_INT>& centered_velocities);
    void Compute_Total_Conserved_Quantity(const bool update_boundary_flux,const T dt,TV_DIMENSION& total_conserved_quantity);
    void Invalidate_Ghost_Cells();
    void Warn_For_Low_Internal_Energy() const;
    bool Equal_Real_Data(const T_ARRAYS_DIMENSION_SCALAR& U1,const T_ARRAYS_DIMENSION_SCALAR& U2) const;
    void Fill_Ghost_Cells(const T dt,const T time,const int ghost_cells) const;
    void Get_Dirichlet_Boundary_Conditions(const T dt,const T time);
    void Advance_One_Time_Step_Forces(const T dt,const T time);
    void Compute_Cavitation_Velocity(ARRAY<T,TV_INT>& rho_n, ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_n, T_ARRAYS_DIMENSION_SCALAR& momentum_n);
    void Advance_One_Time_Step_Explicit_Part(const T dt,const T time,const int rk_substep,const int rk_order);
    void Advance_One_Time_Step_Implicit_Part(const T dt,const T time);
    void Clamp_Internal_Energy(const T dt,const T time); // TODO(kwatra): Do we really need dt, time here?
    void Clamp_Internal_Energy_Ghost(T_ARRAYS_DIMENSION_SCALAR& U_ghost,const int number_of_ghost_cells) const;
    void Remove_Added_Internal_Energy(const T dt,const T time); // TODO(kwatra): Do we really need dt, time here?
    void Euler_Step(const T dt,const T time,const bool is_time_n);
    T CFL_Using_Sound_Speed() const;
    T CFL(const T time) const;
    void Set_Eigensystems(const bool advection_only);
    void Log_Parameters() const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
