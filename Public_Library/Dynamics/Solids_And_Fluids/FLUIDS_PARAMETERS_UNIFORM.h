//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUIDS_PARAMETERS_UNIFORM
//#####################################################################
#ifndef __FLUIDS_PARAMETERS_UNIFORM__    
#define __FLUIDS_PARAMETERS_UNIFORM__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_FORWARD.h>
#include <Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS.h>
namespace PhysBAM{

template<class TV> class MPI_UNIFORM_GRID;
template<class TV> class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM;
template<class TV> class EULER_UNIFORM;
template<class TV> class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES;
template<class TV> class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES;
template<class TV> class LAPLACE_UNIFORM;
template<class TV> class PARTICLE_LEVELSET_EVOLUTION_UNIFORM;

template<class TV>
class FLUIDS_PARAMETERS_UNIFORM:public FLUIDS_PARAMETERS<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T> T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
public:
    typedef FLUIDS_PARAMETERS<TV> BASE;
    using BASE::smoke;using BASE::fire;using BASE::water;using BASE::use_body_force;using BASE::grid;
    using BASE::soot_container;using BASE::soot_fuel_container;using BASE::density_container;using BASE::temperature_container;
    using BASE::use_soot_fuel_combustion;using BASE::burn_temperature_threshold;using BASE::burn_rate;
    using BASE::soot_fuel_calorific_value;
    using BASE::domain_walls;
    using BASE::callbacks;using BASE::gravity;using BASE::phi_boundary;using BASE::fluid_boundary;using BASE::fluid_boundary_water;
    using BASE::phi_boundary_water;using BASE::boundary_mac_slip;using BASE::normal_flame_speed;using BASE::curvature_flame_speed;using BASE::surface_tension;
    using BASE::variable_surface_tension;using BASE::viscosity;using BASE::viscosity_fuel;using BASE::variable_viscosity;using BASE::implicit_viscosity;
    using BASE::implicit_viscosity_iterations;using BASE::use_vorticity_confinement;using BASE::confinement_parameter;
    using BASE::use_vorticity_confinement_fuel;using BASE::use_variable_vorticity_confinement;using BASE::use_strain;using BASE::elastic_modulus;using BASE::plasticity_alpha;
    using BASE::plasticity_gamma;using BASE::adhesion_coefficient;using BASE::confinement_parameter_fuel;using BASE::use_explicit_part_of_implicit_viscosity;
    using BASE::levelset_refinement_bandwidth;using BASE::kolmogorov;using BASE::use_external_velocity;
    using BASE::use_soot;using BASE::use_density;using BASE::use_temperature;using BASE::density;using BASE::soot_advection_order;
    using BASE::density_fuel;using BASE::temperature_products;using BASE::temperature_fuel;using BASE::number_particles_per_cell;using BASE::turbulence_lowest_angular_frequency;
    using BASE::turbulence_update_frame_rate;using BASE::use_non_zero_divergence;using BASE::solve_neumann_regions;using BASE::solve_single_cell_neumann_regions;
    using BASE::temperature_buoyancy_constant;using BASE::object_friction;
    using BASE::move_grid;using BASE::move_grid_explicitly;using BASE::moving_grid_number_of_cells;using BASE::write_velocity;using BASE::write_levelset;using BASE::write_particles;
    using BASE::write_debug_data;using BASE::restart_data_write_rate;using BASE::write_removed_positive_particles;using BASE::write_removed_negative_particles;using BASE::write_strain;
    using BASE::write_ghost_values;using BASE::solid_affects_fluid;using BASE::fluid_affects_solid;using BASE::Initialize_Density_And_Temperature;using BASE::separation_velocity_tolerance;using BASE::write_restart_data;
    using BASE::use_separation_inside_water;using BASE::adhesion_normal_strain;using BASE::modify_wall_tangential_velocities;using BASE::store_particle_ids;using BASE::turbulence;
    using BASE::collision_bodies_affecting_fluid;using BASE::collidable_contour_value;using BASE::collidable_phi_replacement_value;
    using BASE::semi_lagrangian;
    using BASE::phi_boundary_reflection;using BASE::phi_boundary_multiphase;using BASE::flood_fill_for_bubbles;using BASE::adhesion_half_bandwidth;using BASE::sph;
    using BASE::use_maccormack_semi_lagrangian_advection;using BASE::use_maccormack_for_incompressible;using BASE::use_maccormack_for_level_set;using BASE::number_of_ghost_cells;
    using BASE::cfl;using BASE::density_buoyancy_constant;using BASE::rho_top;using BASE::rho_bottom;using BASE::density_buoyancy_threshold;using BASE::use_sph_for_removed_negative_particles;
    using BASE::analytic_test;using BASE::hamilton_jacobi_weno;using BASE::compressible;using BASE::compressible_boundary;using BASE::compressible_pressure_boundary;
    using BASE::compressible_eos;using BASE::compressible_conservation_method;using BASE::compressible_set_max_time_step;using BASE::compressible_max_time_step;
    using BASE::compressible_spatial_order;using BASE::compressible_rungekutta_order;using BASE::compressible_timesplit;using BASE::compressible_apply_isobaric_fix;using BASE::compressible_apply_cavitation_correction;
    using BASE::write_flattened_particles;using BASE::use_poisson;using BASE::simulate;using BASE::use_slip;using BASE::removed_positive_particle_buoyancy_constant;
    using BASE::bandwidth_without_maccormack_near_interface;

    MPI_UNIFORM_GRID<TV>* mpi_grid;
    GRID<TV> p_grid;

    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>* particle_levelset_evolution;
    INCOMPRESSIBLE_UNIFORM<TV>* incompressible;
    PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>* particle_levelset_evolution_multiple;
    INCOMPRESSIBLE_MULTIPHASE_UNIFORM<TV>* incompressible_multiphase;
    SPH_EVOLUTION_UNIFORM<TV>* sph_evolution;
    ARRAY<bool,TV_INT>& maccormack_node_mask;
    ARRAY<bool,TV_INT>& maccormack_cell_mask;
    ARRAY<bool,FACE_INDEX<TV::m> >& maccormack_face_mask;
    ADVECTION_MACCORMACK_UNIFORM<TV,T,T_ADVECTION_SEMI_LAGRANGIAN_SCALAR>& maccormack_semi_lagrangian;
    EULER_UNIFORM<TV>* euler;
    SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>* euler_solid_fluid_coupling_utilities;
    COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>* compressible_incompressible_coupling_utilities;
    PROJECTION_DYNAMICS_UNIFORM<TV>* projection;

    // multiphase parameters
    ARRAY<T> masses; // used to keep track of mass loss in the driver
    ARRAY<T> densities;
    ARRAY<T> viscosities;
    ARRAY<T,VECTOR<int,2> > surface_tensions;
    ARRAY<T> confinement_parameters;
    ARRAY<bool> dirichlet_regions;
    ARRAY<bool> pseudo_dirichlet_regions;
    ARRAY<bool> fuel_region; // is this region a fuel region
    bool use_reacting_flow;
    ARRAY<T,VECTOR<int,2> > normal_flame_speeds; // specify 0 for non reacting flow.  must be symmetric
    ARRAY<T,VECTOR<int,2> > curvature_flame_speeds; // specify 0 for non reacting flow.  must be skew symmetric
    int number_of_regions;
    bool use_flame_speed_multiplier;
    bool use_dsd;
    ARRAY<bool> use_multiphase_strain;
    ARRAY<T> elastic_moduli;
    ARRAY<T> plasticity_alphas;
    ARRAY<T> plasticity_gammas;
    bool use_psi_R;
    bool use_levelset_viscosity;
    bool print_viscosity_matrix;
    bool use_second_order_pressure;
    bool use_surface_solve;
    int projection_scale;
    VECTOR<bool,TV::m> periodic_boundary;

    FLUIDS_PARAMETERS_UNIFORM(const int number_of_regions,const typename FLUIDS_PARAMETERS<TV>::TYPE type);
    virtual ~FLUIDS_PARAMETERS_UNIFORM();

//#####################################################################
    virtual void Initialize_Grids();
    void Initialize_Fluid_Evolution(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void Use_Fluid_Coupling_Defaults() override;
    void Use_No_Fluid_Coupling_Defaults() override;
    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
    void Delete_Particles_Inside_Objects(const T time);
    void Initialize_Number_Of_Regions(const int number_of_regions_input);
private:
    template<class T_PARTICLES> void Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*,TV_INT>& particles,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
public:
    void Set_Projection(PROJECTION_DYNAMICS_UNIFORM<TV>* projection_input);
    void Update_Fluid_Parameters(const T dt,const T time);
    void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time);
    void Apply_Isobaric_Fix(const T dt,const T time);
    void Get_Neumann_And_Dirichlet_Boundary_Conditions(LAPLACE_UNIFORM<TV>* elliptic_solver,
            ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    void Set_Domain_Boundary_Conditions(LAPLACE_UNIFORM<TV>& elliptic_solver,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time);
    void Blend_In_External_Velocity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    void Move_Grid(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time);
    void Move_Grid(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const TV_INT& shift_domain,const T time);
    void Adjust_Strain_For_Object(LEVELSET<TV>& levelset_object,ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& e_ghost,const T time);
    void Combustion(const T dt,const T time);
    void Evolve_Soot(const T dt,const T time);
    template<class T_ARRAYS_PARTICLES> int Total_Number_Of_Particles(const T_ARRAYS_PARTICLES& particles) const;
    template<class T_ARRAYS_PARTICLES> void Write_Particles(const STREAM_TYPE stream_type,const PARTICLES<TV>& template_particles,const T_ARRAYS_PARTICLES& particles,
        const std::string& output_directory,const std::string& prefix,const int frame) const;
    template<class T_PARTICLES,class T_ARRAYS_PARTICLES> void Read_Particles(const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,
        const std::string& output_directory,const std::string& prefix,const int frame);
    void Read_Output_Files(const std::string& output_directory,const int frame);
    void Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int first_frame,const int frame) const;
    void Log_Parameters() const override;
//#####################################################################
};      
}
#endif
