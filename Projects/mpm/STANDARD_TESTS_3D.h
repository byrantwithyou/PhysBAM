//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_3D__
#define __STANDARD_TESTS_3D__

#include <Tools/Parsing/PARSE_ARGS.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include "STANDARD_TESTS_BASE.h"

namespace PhysBAM{

template<class TV> class STANDARD_TESTS;
template<class T> class TRIANGULATED_SURFACE;
template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class T> class CYLINDER;

template<class T>
class STANDARD_TESTS<VECTOR<T,3> >:public STANDARD_TESTS_BASE<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef STANDARD_TESTS_BASE<TV> BASE;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;

public:
    using BASE::last_frame;using BASE::grid;using BASE::particles;
    using BASE::frame_title;using BASE::write_substeps_level;using BASE::particles_per_cell;
    using BASE::substeps_delay_frame;using BASE::scale_mass;using BASE::scale_E;
    using BASE::viewer_dir;using BASE::restart;using BASE::dt;using BASE::time;
    using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;using BASE::order;
    using BASE::ghost;using BASE::use_affine;
    using BASE::use_midpoint;using BASE::print_stats;using BASE::flip;using BASE::cfl;
    using BASE::newton_tolerance;using BASE::newton_iterations;using BASE::solver_tolerance;
    using BASE::lagrangian_forces;using BASE::mass;using BASE::gather_scatter;
    using BASE::solver_iterations;using BASE::test_diff;using BASE::threads;using BASE::test_number;
    using BASE::resolution;using BASE::Seed_Particles;using BASE::Add_Gravity;
    using BASE::Add_Fixed_Corotated;using BASE::random;using BASE::Seed_Lagrangian_Particles;using BASE::Add_Clamped_Plasticity;
    using BASE::Add_Force;using BASE::Add_Walls;using BASE::data_directory;
    using BASE::stream_type;using BASE::use_oldroyd;using BASE::force_helper;
    using BASE::Add_Neo_Hookean;using BASE::Add_St_Venant_Kirchhoff_Hencky_Strain;using BASE::Add_Collision_Object;
    using BASE::Add_Particle;using BASE::Add_Lambda_Particles;using BASE::Add_Penalty_Collision_Object;
    using BASE::penalty_collisions_stiffness;using BASE::penalty_damping_stiffness;
    using BASE::quad_F_coeff;using BASE::use_penalty_collisions;using BASE::Add_Drucker_Prager;
    using BASE::use_implicit_plasticity;using BASE::no_implicit_plasticity;
    using BASE::use_theta_c;using BASE::theta_c;using BASE::use_theta_s;using BASE::theta_s;
    using BASE::hardening_factor;using BASE::use_hardening_factor;using BASE::max_hardening;using BASE::use_max_hardening;
    using BASE::hardening_mast_case;using BASE::use_hardening_mast_case;
    using BASE::m;using BASE::s;using BASE::kg;using BASE::unit_p;using BASE::unit_mu;using BASE::unit_rho;
    using BASE::forced_collision_type;using BASE::scale_speed;using BASE::Add_Drucker_Prager_Case;
    using BASE::friction;using BASE::friction_is_set;using BASE::sigma_Y;using BASE::use_cohesion;using BASE::destroy;
    using BASE::collision_objects;using BASE::Set_Lame_On_Particles;using BASE::plasticity_models;
    using BASE::write_output_files;using BASE::read_output_files;
    using BASE::dump_collision_objects;using BASE::Perturb;using BASE::Uniform;using BASE::extra_T;using BASE::extra_int;
    using BASE::Set_Grid;using BASE::Add_Callbacks;
    using BASE::sph_rel;using BASE::Add_Source;
    

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS();

    void Initialize() override;

    // additional storage
    int foo_int1;
    T foo_T1;
    T foo_T2;
    T foo_T3;
    T foo_T4;
    T foo_T5;
    bool use_foo_T1;
    bool use_foo_T2;
    bool use_foo_T3;
    bool use_foo_T4;
    bool use_foo_T5;

//#####################################################################
};
}

#endif
