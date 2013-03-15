//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURFACE_TENSION
//#####################################################################
//   1. Circle with surface tension, pressure jump condition
//   2. Surface tension on a thin solid
//   3. Circle with surface tension, pressure jump condition, coupled to passive solid
//   4. Oscillating deformed circle
//   5. Oscillating deformed circle with massless solid
//   6. Analytic viscosity test
//   7. Two-phase rising bubble test
//   8. Analytic curvature test
//   9. Stationary circle test for JCP with pressure jump version
//   10. Stationary circle test for JCP with phi and front tracking version
//   11. Small rising bubble test
//   12. Popinet example
//   13. Small rising bubble test using pressure jump
//#####################################################################
#ifndef __SURFACE_TENSION__
#define __SURFACE_TENSION__

#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_POINT_ATTRACTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_CUT.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_PHI.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
#include <PhysBAM_Dynamics/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <fstream>
namespace PhysBAM{
template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);
template<class TV> void Add_Debug_Particle(const TV& X){Add_Debug_Particle(X,VECTOR<typename TV::SCALAR,3>(1,0,0));}
template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);

template<class T_input>
class SURFACE_TENSION:public PLS_FSI_EXAMPLE<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef ARRAY<T,FACE_INDEX<2> > T_FACE_ARRAYS_SCALAR;
    typedef ARRAY<bool,FACE_INDEX<2> > T_FACE_ARRAYS_BOOL;
public:
    typedef PLS_FSI_EXAMPLE<TV> BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::solids_fluids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::two_phase;
    using BASE::parse_args;using BASE::test_number;using BASE::resolution;using BASE::data_directory;using BASE::convection_order;using BASE::use_pls_evolution_for_structure;
    using BASE::Mark_Outside;using BASE::Add_Thin_Shell_To_Fluid_Simulation;using BASE::Add_To_Fluid_Simulation;

    SOLIDS_STANDARD_TESTS<TV> solids_tests;

    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_matrix;
    bool print_each_matrix;
    bool use_full_ic;
    bool output_iterators;
    bool use_decoupled_viscosity;
    T max_dt;
    T exact_dt;
    T current_dt;
    bool implicit_solid,no_implicit_solid,use_cut_volume,use_low_order_advection;

    GEOMETRY_PARTICLES<TV> debug_particles;
    SEGMENTED_CURVE_2D<T>* front_tracked_structure;
    SEGMENTED_CURVE_2D<T>* rebuild_curve;
    ARRAY<TV> saved_tracked_particles_X;
    DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* deformable_collisions;
    FLUID_TO_SOLID_INTERPOLATION_CUT<TV>* fsi;

    int number_surface_particles;
    bool rebuild_surface;
    ARRAY<VECTOR<int,2> > particle_segments;
    FREE_PARTICLES<TV>* free_particles;
    ARRAY<typename MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED<TV>::ENTRY> fluid_interpolation_entries;
    ARRAY<int> solid_interpolation_entries;
    ARRAY<bool,TV_INT>* psi_D;

    T circle_radius;
    T circle_perturbation;
    int oscillation_mode;
    bool use_massless_structure;
    ARRAY<int>* coupled_particles;
    bool make_ellipse,use_phi,remesh;
    T m,s,kg;
    int solid_refinement;
    T solid_density,solid_width,analytic_solution;
    T omega;
    T laplace_number,surface_tension;
    bool use_T_nu;
    bool use_viscosity;
    T linear_force,rand;
    int cut_order;

    SURFACE_TENSION(const STREAM_TYPE stream_type);
    virtual ~SURFACE_TENSION();

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    //void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Reseed_Mask(ARRAY<bool,TV_INT>*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE {return false;}
    void Mark_Outside(ARRAY<bool,TV_INT>& outside) {}

//#####################################################################
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE;
    void Register_Options() PHYSBAM_OVERRIDE;
    void Parse_Options() PHYSBAM_OVERRIDE;
    void Add_Rigid_Body_Walls(const T coefficient_of_restitution=(T).5,const T coefficient_of_friction=(T).5,ARRAY<int>* walls_added=0);
    void Parse_Late_Options() PHYSBAM_OVERRIDE;
    void Initialize_Advection() PHYSBAM_OVERRIDE;
    void Initialize_Phi() PHYSBAM_OVERRIDE;
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE;
    void Initialize_Velocities() PHYSBAM_OVERRIDE;
    void Set_Dirichlet_Boundary_Conditions(const T time);
    void Get_Source_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_BOOL& psi_N,const T time) PHYSBAM_OVERRIDE;
    void Initialize_Bodies() PHYSBAM_OVERRIDE;
    void Kang_Circle(bool use_surface);
    void Oscillating_Circle(bool use_surface);
    void Test_Analytic_Velocity(T time);
    void Test_Analytic_Pressure(T time);
    void Solid_Circle();
    void Adjust_Phi_With_Objects(const T time);
    void Sync_Particle_To_Level_Set(int p);
    void Sync_Front_Tracked_Particles_To_Level_Set();
    void Divide_Segment(int e);
    void Swap_Particles(int p,int r);
    void Swap_Segments(int e,int f);
    void Remove_Particle(int p);
    T Compute_New_Mass(int p);
    void Copy_Front_Tracked_Velocity_From_Fluid();
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE;
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE;
    void Write_Output_Files(const int frame) const;
    void Initialize_Surface_Particles(int number);
    void Rebuild_Surface();
    void Substitute_Coupling_Matrices(KRYLOV_SYSTEM_BASE<T>& coupled_system,T dt,T current_velocity_time,T current_position_time,bool velocity_update,bool leakproof_solve) PHYSBAM_OVERRIDE;
    void Advance_One_Time_Step_Begin_Callback(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE;
    static GEOMETRY_PARTICLES<TV>*  Store_Debug_Particles(GEOMETRY_PARTICLES<TV>* particle=0);
    void FSI_Analytic_Test();
    void Sine_Wave();
    void Initialize_Sine_Phi();
//#####################################################################
};

}
#endif

