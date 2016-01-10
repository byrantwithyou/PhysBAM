//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_BASE__
#define __STANDARD_TESTS_BASE__

#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Deformables/Standard_Tests/DEFORMABLES_STANDARD_TESTS.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <functional>

namespace PhysBAM{

template<class TV> class STANDARD_TESTS;
template<class T,int d> class ISOTROPIC_CONSTITUTIVE_MODEL;

template<class TV>
class STANDARD_TESTS_BASE:public MPM_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef MPM_EXAMPLE<TV> BASE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_VOLUME;

public:
    using BASE::initial_time;using BASE::last_frame;using BASE::grid;using BASE::particles;
    using BASE::mass;using BASE::force_helper;using BASE::use_oldroyd;
    using BASE::frame_title;using BASE::write_substeps_level;using BASE::gather_scatter;
    using BASE::collision_objects;using BASE::substeps_delay_frame;
    using BASE::output_directory;using BASE::mass_contour;using BASE::plasticity_models;
    using BASE::restart;using BASE::dt;using BASE::time;using BASE::use_early_gradient_transfer;
    using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;
    using BASE::ghost;using BASE::use_reduced_rasterization;using BASE::use_affine;;using BASE::use_f2p;
    using BASE::use_midpoint;using BASE::use_symplectic_euler;
    using BASE::print_stats;using BASE::only_write_particles;
    using BASE::flip;using BASE::cfl;using BASE::newton_tolerance;
    using BASE::newton_iterations;using BASE::solver_tolerance;using BASE::solver_iterations;
    using BASE::test_diff;using BASE::threads;using BASE::weights;using BASE::incompressible;
    using BASE::lagrangian_forces;using BASE::stream_type;
    using BASE::Add_Force;using BASE::Set_Weights;using BASE::deformable_body_collection;
    using BASE::Add_Collision_Object;using typename BASE::COLLISION_TYPE;using BASE::data_directory;
    using BASE::kkt;using BASE::Add_Fluid_Wall;using BASE::quad_F_coeff;

    int test_number;
    int resolution;
    bool user_resolution;
    int stored_last_frame;
    bool user_last_frame;
    int order;
    int seed;
    int particles_per_cell;
    bool regular_seeding,no_regular_seeding;
    T scale_mass;
    T scale_E;
    T scale_speed;
    T penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length;
    T penalty_damping_stiffness;
    bool use_penalty_collisions;
    bool use_plasticity,use_theta_c,use_theta_s,use_hardening_factor,use_max_hardening;
    T theta_c,theta_s,hardening_factor,max_hardening;
    bool use_implicit_plasticity,no_implicit_plasticity;
    int hardening_mast_case;
    bool use_hardening_mast_case;
    bool override_output_directory;
    T m,s,kg;
    T unit_p,unit_rho,unit_mu;
    int forced_collision_type;
    T friction;
    bool friction_is_set;
    T sigma_Y;
    bool use_cohesion;
    std::function<void (int frame)> write_output_files;
    std::function<void (int frame)> read_output_files;
    std::function<void (int frame)> begin_frame;
    std::function<void (int frame)> end_frame;
    std::function<void (T time)> begin_time_step;
    std::function<void (T time)> end_time_step;
    std::function<void ()> destroy;

    RANDOM_NUMBERS<T> random;
    DEFORMABLES_STANDARD_TESTS<TV> tests;

    STANDARD_TESTS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_BASE();

    void Seed_Particles_Poisson(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,int particles_per_cell);

    template<class T_OBJECT> typename enable_if<!is_base_of<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>::type
    Seed_Particles_Poisson(const T_OBJECT& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,int particles_per_cell)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles_Poisson(obj,V,dV,density,particles_per_cell);}

    void Seed_Particles_Uniform(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,const GRID<TV>& seed_grid);

    template<class T_OBJECT> typename enable_if<!is_base_of<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>::type
    Seed_Particles_Uniform(const T_OBJECT& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,const GRID<TV>& seed_grid)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles_Uniform(obj,V,dV,density,seed_grid);}

    void Seed_Particles(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,int particles_per_cell);

    template<class T_OBJECT> typename enable_if<!is_base_of<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>::type
    Seed_Particles(const T_OBJECT& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,int particles_per_cell)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles(obj,V,dV,density,particles_per_cell);}

    template<class T_STRUCTURE>
    T_STRUCTURE& Seed_Lagrangian_Particles(T_STRUCTURE& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,bool use_constant_mass,bool destroy_after=true);

    void Add_Penalty_Collision_Object(IMPLICIT_OBJECT<TV>* io,const T coefficient_of_friction=(T)0);
    template<class OBJECT> typename enable_if<!is_pointer<OBJECT>::value>::type
    Add_Penalty_Collision_Object(const OBJECT& object,const T coefficient_of_friction=(T)0)
    {Add_Penalty_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<OBJECT>(object),coefficient_of_friction);}

    void Add_Particle(const TV& X,std::function<TV(const TV&)> V,std::function<MATRIX<T,TV::m>(const TV&)> dV,
        const T mass,const T volume);
    void Add_Lambda_Particles(ARRAY<int>* affected_particles,T E,T nu,T density,bool no_mu=false,T porosity=(T)1,T saturation_level=(T)1);
    int Add_Gravity(TV g,ARRAY<int>* affected_particles=0);
    int Add_Fixed_Corotated(T E,T nu,ARRAY<int>* affected_particles=0,bool no_mu=false);
    int Add_Neo_Hookean(T E,T nu,ARRAY<int>* affected_particles=0);
    int Add_Fixed_Corotated(T_VOLUME& object,T E,T nu);
    int Add_Neo_Hookean(T_VOLUME& object,T E,T nu);
    int Add_St_Venant_Kirchhoff_Hencky_Strain(T E,T nu,ARRAY<int>* affected_particles=0,bool no_mu=false);
    int Add_Drucker_Prager(T E,T nu,T a0,T a1,T a3,T a4,ARRAY<int>* affected_particles=0,bool no_mu=false,T sigma_Y=0);
    int Add_Drucker_Prager(T E,T nu,T phi_F,ARRAY<int>* affected_particles=0,bool no_mu=false,T sigma_Y=0);
    int Add_Drucker_Prager_Case(T E,T nu,int case_num_F,ARRAY<int>* affected_particles=0,bool no_mu=false);
    void Add_Walls(int flags,COLLISION_TYPE type,T friction,T inset,bool penalty); // -x +x -y +y [ -z +z ], as bit flags
    int Add_Clamped_Plasticity(ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& icm,T theta_c,T theta_s,T max_hardening,
        T hardening_factor,ARRAY<int>* affected_particles);
    void Set_Lame_On_Particles(T E,T nu,ARRAY<int>* affected_particles=0);
    void Update_Variable_Lame_Parameters_On_Constitutive_Models();
//#####################################################################
};
}

#endif
