//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KANG
//#####################################################################
//   1. Circle with surface tension, pressure jump condition
//   2. Oscillating deformed circle
//   3. Two-phase rising bubble test
//   4. Poisson Test
//   5. Poiseuille Flow Test
//   6. Couette Flow Test
//   7. Circular Couette Flow Test
//   8. Radial Flow Test
//   9. Kang Example 1
//#####################################################################
#ifndef __KANG__
#define __KANG__

#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
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
#include <Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION.h>
#include <Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_CUT.h>
#include <Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_PHI.h>
#include <Dynamics/Coupled_Evolution/MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED.h>
#include <Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <Dynamics/Coupled_Evolution/SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
#include <Dynamics/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <fstream>
namespace PhysBAM{
template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);
template<class TV> void Add_Debug_Particle(const TV& X){Add_Debug_Particle(X,VECTOR<typename TV::SCALAR,3>(1,0,0));}
template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);

template<class T_input>
class KANG:public PLS_FSI_EXAMPLE<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    typedef PLS_FSI_EXAMPLE<TV> BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::solids_fluids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Add_Volumetric_Body_To_Fluid_Simulation;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::two_phase;
    using BASE::test_number;using BASE::resolution;using BASE::data_directory;using BASE::convection_order;using BASE::use_pls_evolution_for_structure;
    using BASE::Mark_Outside;using BASE::use_kang;using BASE::print_matrix;using BASE::test_system;
    using BASE::m;using BASE::s;using BASE::kg;

    SOLIDS_STANDARD_TESTS<TV> solids_tests;

    bool output_iterators;
    T max_dt;
    T exact_dt;
    bool implicit_solid;

    GEOMETRY_PARTICLES<TV> debug_particles;

    T circle_radius;
    T circle_perturbation;
    int oscillation_mode;
    bool make_ellipse;
    T omega;
    T laplace_number,surface_tension;

    T uleft,uright;

    // for tests 7 and 8
    T r_n,r_I,r_p,u_n0,u_p0;

    KANG(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~KANG();

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
    void After_Initialization() PHYSBAM_OVERRIDE;
    void Initialize_Advection() PHYSBAM_OVERRIDE;
    void Initialize_Phi() PHYSBAM_OVERRIDE;
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE;
    void Initialize_Velocities() PHYSBAM_OVERRIDE;
    void Set_Dirichlet_Boundary_Conditions(const T time);
    void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) PHYSBAM_OVERRIDE;
    void Initialize_Bodies() PHYSBAM_OVERRIDE;
    void Kang_Circle();
    void Poisson_Test();
    void Poiseuille_Flow_Test();
    void Couette_Flow_Test();
    void Circular_Couette_Flow_Test();
    void Radial_Flow_Test();
    void Kang_Example_1();
    void Oscillating_Circle();
    void Test_Analytic_Velocity(T time);
    void Test_Analytic_Pressure(T time);
    //void Solid_Circle();
    //void Sync_Particle_To_Level_Set(int p);
    //void Sync_Front_Tracked_Particles_To_Level_Set();
    //void Divide_Segment(int e);
    //void Swap_Particles(int p,int r);
    //void Swap_Segments(int e,int f);
    //void Remove_Particle(int p);
    //T Compute_New_Mass(int p);
    //void Copy_Front_Tracked_Velocity_From_Fluid();
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE;
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE;
    void Write_Output_Files(const int frame) const;
    //void Initialize_Surface_Particles(int number);
    //void Rebuild_Surface();
    void Advance_One_Time_Step_Begin_Callback(const T dt,const T time) PHYSBAM_OVERRIDE;
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE;
    static GEOMETRY_PARTICLES<TV>*  Store_Debug_Particles(GEOMETRY_PARTICLES<TV>* particle=0);
    //void FSI_Analytic_Test();
    //void Sine_Wave();
    //void Initialize_Sine_Phi();
    void Set_Analytic_Velocity(const T time,ARRAY<T,FACE_INDEX<TV::dimension> >& u) const;
    void Set_Boundary_Conditions_Callback(ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& psi_D_value,
        ARRAY<T,FACE_INDEX<TV::dimension> >& psi_N_value) const PHYSBAM_OVERRIDE;
//#####################################################################
};

}
#endif

