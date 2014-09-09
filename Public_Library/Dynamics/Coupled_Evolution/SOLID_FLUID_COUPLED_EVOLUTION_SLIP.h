//#####################################################################
// Copyright 2008-2009, Elliot English, Nipun Kwatra, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_COUPLED_EVOLUTION_SLIP
//#####################################################################
#ifndef __SOLID_FLUID_COUPLED_EVOLUTION_SLIP__
#define __SOLID_FLUID_COUPLED_EVOLUTION_SLIP__

#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS_BINARY_UNIFORM.h>
#include <Geometry/Topology/TOPOLOGY_POLICY.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <Dynamics/Coupled_Evolution/COUPLED_SYSTEM_VECTOR.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
namespace PhysBAM{

template<class T> class FRACTURE_PATTERN;
template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class FLUIDS_PARAMETERS_UNIFORM;
template<class TV> class POISSON_COLLIDABLE_UNIFORM;
template<class TV> class GENERALIZED_VELOCITY;
template<class TV> class GENERALIZED_MASS;
template<class TV> class MATRIX_SOLID_INTERPOLATION;
template<class TV> class SOLIDS_FLUIDS_PARAMETERS;
template<class TV> class BACKWARD_EULER_SYSTEM;
template<class TV> class FLUID_COLLECTION;
template<class TV> class EULER_PROJECTION_UNIFORM;
template<class TV> class IMPLICIT_BOUNDARY_CONDITION_COLLECTION;
template<class TV> class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO;

template<class TV>
class SOLID_FLUID_COUPLED_EVOLUTION_SLIP:public NEWMARK_EVOLUTION<TV>,public PROJECTION_DYNAMICS_UNIFORM<TV>
{
public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
    typedef ARRAY<int,FACE_INDEX<TV::m> > T_FACE_ARRAYS_INT;

protected:
    typedef NEWMARK_EVOLUTION<TV> BASE;
    typedef PROJECTION_DYNAMICS_UNIFORM<TV> FLUID_BASE;
    using BASE::solid_body_collection;using BASE::solids_parameters;using BASE::B_full;using BASE::rigid_B_full;
    using BASE::repulsions;using BASE::rigids_evolution_callbacks;using BASE::rigid_body_collisions;
    using FLUID_BASE::p;using FLUID_BASE::poisson;using BASE::Prepare_Backward_Euler_System;using BASE::krylov_vectors;

    COUPLED_SYSTEM_VECTOR<TV> coupled_x,coupled_b;

    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies;
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters;
    FLUID_COLLECTION<TV>& fluid_collection;
    ARRAY<TV> leakproof_empty_V,temp_solid_full,second_temp_solid_full;
    ARRAY<TWIST<TV> > rigid_temp_solid_full,rigid_second_temp_solid_full,leakproof_empty_rigid_V;
    
    T divergence_scaling;
    bool disable_thinshell;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> coupled_vectors;

public:
    UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& iterator_info;
    IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>& boundary_condition_collection;
protected:

    GRID<TV>& grid;
    ARRAY<T,FACE_INDEX<TV::m> > fluids_face_velocities; // Stores combined compressible-incompressible face velocities.
public:
    ARRAY<T,TV_INT> pressure;
    ARRAY<T,TV_INT> density;
    ARRAY<bool,FACE_INDEX<TV::m> > solved_faces;
protected:
    ARRAY<TV,TV_INT> centered_velocity;
    ARRAY<T,TV_INT> one_over_rho_c_squared;
    ARRAY<T,TV_INT> p_advected_over_rho_c_squared_dt;
    ARRAY<T,TV_INT> p_advected;

public:
    FRACTURE_PATTERN<T>* fracture_pattern;
private:
    ARRAY<TV> pressure_impulses;
    ARRAY<TWIST<TV> > pressure_impulses_twist;

    ARRAY<bool,FACE_INDEX<TV::m> > cached_psi_N;
    ARRAY<T,COUPLING_CONSTRAINT_ID> coupling_face_velocities_cached;
    COUPLING_CONSTRAINT_ID number_of_coupling_faces_cached;
    ARRAY<FACE_INDEX<TV::dimension> > indexed_faces_cached;
    T time_cached;
    bool cached_coupling_face_data;
public:
    using BASE::print_matrix;
    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_rhs;
    bool print_each_matrix;
    bool output_iterators;
    bool use_viscous_forces;
    bool two_phase;
    bool use_full_ic;

    SOLID_FLUID_COUPLED_EVOLUTION_SLIP(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,
        EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities_input,FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_input,
        SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters_input,FLUID_COLLECTION<TV>& fluid_collection_input);
    virtual ~SOLID_FLUID_COUPLED_EVOLUTION_SLIP();

    bool Cell_To_Cell_Visible(const int axis,const TV_INT& first_cell,const TV_INT& second_cell) const
    {assert(first_cell==second_cell+TV_INT::Axis_Vector(axis));return collision_bodies.cell_neighbors_visible(second_cell)(axis);}

//#####################################################################
    void Initialize_Grid_Arrays();
    bool Simulate_Incompressible_Fluids() const;
    bool Simulate_Compressible_Fluids() const;
    bool Simulate_Fluids() const;
    bool Simulate_Solids() const;
    void Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update) PHYSBAM_OVERRIDE;
    void Process_Collisions(const T dt,const T time,const bool advance_rigid_bodies) PHYSBAM_OVERRIDE;
    void Apply_Pressure(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,bool scale_by_dt=false);
    BACKWARD_EULER_SYSTEM<TV>* Setup_Solids(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update,const bool leakproof_solve);
    void Setup_Fluids(ARRAY<T,FACE_INDEX<TV::m> >& incompressible_face_velocities,const T current_position_time,const T dt,const bool leakproof_solve);
    void Solve(ARRAY<T,FACE_INDEX<TV::m> >& incompressible_face_velocities,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update,const bool leakproof_solve);
    void Make_Divergence_Free(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    void Apply_Second_Order_Cut_Cell_Method(const T_ARRAYS_INT& cell_index_to_divergence_matrix_index,const T_FACE_ARRAYS_INT& face_index_to_matrix_index,ARRAY<T>& b);
    void Apply_Viscosity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time);
    void Setup_Boundary_Condition_Collection();
private:
    void Warn_For_Exposed_Dirichlet_Cell(const ARRAY<bool,TV_INT>& psi_D,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N);
    void Set_Cached_Psi_N_And_Coupled_Face_Data(const COLLISION_AWARE_INDEX_MAP<TV>& index_map,
        const MATRIX_SOLID_INTERPOLATION<TV>& solid_interpolation,const T time);
    void Fill_Coupled_Face_Data(const COUPLING_CONSTRAINT_ID number_of_coupling_faces,const ARRAY<FACE_INDEX<TV::dimension> >& indexed_faces,
        const ARRAY<T,COUPLING_CONSTRAINT_ID>& coupling_face_data,ARRAY<T,FACE_INDEX<TV::m> >& face_data);
    void Get_Coupled_Faces_And_Interpolated_Solid_Velocities(const COLLISION_AWARE_INDEX_MAP<TV>& index_map,
        const MATRIX_SOLID_INTERPOLATION<TV>& solid_interpolation,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N_domain,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
        ARRAY<T,COUPLING_CONSTRAINT_ID>& coupling_face_velocities);
public:
    void Get_Coupled_Faces_And_Interpolated_Solid_Velocities(const T time,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void Output_Iterators(const STREAM_TYPE stream_type,const char* output_directory,int frame) const;
//#####################################################################
};
}
#endif
