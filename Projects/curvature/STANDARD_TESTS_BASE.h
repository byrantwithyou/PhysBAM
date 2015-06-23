//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_BASE
//#####################################################################
#ifndef __STANDARD_TESTS_BASE__
#define __STANDARD_TESTS_BASE__

#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <fstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{

template<class TV>
class STANDARD_TESTS_BASE:public SOLIDS_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_OBJECT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::solid_body_collection;
    using BASE::stream_type;using BASE::solids_evolution;using BASE::test_number;using BASE::data_directory;using BASE::m;using BASE::s;using BASE::kg;

    SOLIDS_STANDARD_TESTS<TV> tests;

    bool test_forces;
    ARRAY<INTERPOLATION_CURVE<T,FRAME<TV> > > curves;
    ARRAY<int> kinematic_points;
    ARRAY<INTERPOLATION_CURVE<T,TV> > point_curves;
    bool print_matrix;
    int resolution;
    T stiffness_multiplier;
    T thickness_multiplier;
    T curvature_stiffness_multiplier;
    T damping_multiplier;
    ARRAY<int> externally_forced;
    ARRAY<int> constrained_particles;
    ARRAY<TV> constrained_velocities;
    T input_poissons_ratio,input_youngs_modulus;
    T input_friction;
    T ether_drag;
    int rand_seed;
    bool use_rand_seed;
    RANDOM_NUMBERS<T> rand;
    bool use_newmark,use_newmark_be;
    bool project_nullspace;
    BACKWARD_EULER_EVOLUTION<TV>* backward_euler_evolution;
    bool use_penalty_collisions;
    bool use_constraint_collisions;
    bool no_line_search;
    bool no_descent;
    T penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length;
    bool enforce_definiteness;
    T unit_rho,unit_p,unit_N,unit_J;
    T density;
    bool use_penalty_self_collisions;
    T save_dt;
    bool self_collide_surface_only;
    bool use_vanilla_newton;
    int gauss_order;
    int threads;

    STANDARD_TESTS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_BASE();

    void Set_Number_Of_Threads(int threads);
    void Post_Initialization() override;
    void Postprocess_Solids_Substep(const T time,const int substep) override;
    void Apply_Constraints(const T dt,const T time) override;
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) override;
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) override;
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) override;
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) override;
    void Limit_Solids_Dt(T& dt,const T time) override;
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override;
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) override;
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override;
    void Align_Deformable_Bodies_With_Rigid_Bodies() override;
    void Preprocess_Solids_Substep(const T time,const int substep) override;
    void Update_Solids_Parameters(const T time) override;
    void Self_Collisions_Begin_Callback(const T time,const int substep) override;
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) override;
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) override;
    void Preprocess_Frame(const int frame) override;
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override;
    void Update_Time_Varying_Material_Properties(const T time) override;
    void Postprocess_Substep(const T dt,const T time) override;
    void Postprocess_Frame(const int frame) override;
    void Get_Initial_Data_After(bool automatically_add_to_collision_structures);
    void Initialize_Bodies_After();
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override;
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override;
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override;
    void Read_Output_Files_Solids(const int frame) override;
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id);
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id);
    void Preprocess_Substep(const T dt,const T time) override;
    void Add_Constitutive_Model(T_OBJECT& object,T stiffness,T poissons_ratio,T damping);
    GRAVITY<TV>& Add_Gravity();
};
}
#endif
