//#####################################################################
// Copyright.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Circle with surface tension, pressure jump condition
//   2. Surface tension on a thin solid
//   3. Circle with surface tension, pressure jump condition, coupled to passive solid
//   4. Oscillating deformed circle
//   5. Oscillating deformed circle with massless solid
//   6. Analytic viscosity test
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_2D.h>
#include <Deformables/Forces/LINEAR_POINT_ATTRACTION.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <Deformables/Particles/FREE_PARTICLES.h>
#include <Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION.h>
#include <Dynamics/Coupled_Evolution/MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED.h>
#include <Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <Dynamics/Coupled_Evolution/SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
#include <Dynamics/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <fstream>
namespace PhysBAM{
template <class TV> class FLUID_TO_SOLID_INTERPOLATION_CUT;
template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    typedef ARRAY<T,FACE_INDEX<2> > ARRAY<T,FACE_INDEX<TV::m> >;
    typedef ARRAY<bool,FACE_INDEX<2> > ARRAY<bool,FACE_INDEX<TV::m> >;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::solids_fluids_parameters;using BASE::viewer_dir;using BASE::last_frame;using BASE::frame_rate;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::solid_body_collection;using BASE::solids_evolution;
    using BASE::test_number;using BASE::resolution;using BASE::data_directory;

    SOLIDS_STANDARD_TESTS<TV> solids_tests;

    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_matrix;
    bool print_each_matrix;
    bool use_full_ic;
    bool output_iterators;
    bool use_viscous_forces;
    T max_dt;
    T exact_dt;
    T current_dt;
    bool implicit_solid;

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
    bool make_ellipse;
    T m,s,kg;
    int solid_refinement;
    T solid_density,solid_width,analytic_solution,linear_force,rand;
    bool use_viscosity;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS();

//#####################################################################
    void Postprocess_Substep(const T dt,const T time) override;
    void Postprocess_Frame(const int frame) override;
    void Initialize_Advection() override;
    void Initialize_Phi() override;
    void Preprocess_Substep(const T dt,const T time) override;
    void Preprocess_Frame(const int frame) override;
    void Initialize_Velocities() override;
    void Set_Dirichlet_Boundary_Conditions(const T time);
    void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) override;
    void Initialize_Bodies() override;
    void Kang_Circle(bool use_surface);
    void Oscillating_Circle(bool use_surface);
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
    void Limit_Dt(T& dt,const T time) override;
    void Limit_Solids_Dt(T& dt,const T time) override;
    void Write_Output_Files() const;
    void Initialize_Surface_Particles(int number);
    void Rebuild_Surface();
    void Substitute_Coupling_Matrices(KRYLOV_SYSTEM_BASE<T>& coupled_system,T dt,T current_velocity_time,T current_position_time,bool velocity_update,bool leakproof_solve) override;
    void Update_Time_Varying_Material_Properties(const T time) override;
    void FSI_Analytic_Test();
//#####################################################################
};
}
#endif

