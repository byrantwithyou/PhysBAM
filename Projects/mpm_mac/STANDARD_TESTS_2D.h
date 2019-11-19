//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_2D__
#define __STANDARD_TESTS_2D__

#include <Tools/Parsing/PARSE_ARGS.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include "STANDARD_TESTS_BASE.h"

namespace PhysBAM{

template<class TV> class STANDARD_TESTS;

template<class T>
class STANDARD_TESTS<VECTOR<T,2> >:public STANDARD_TESTS_BASE<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef STANDARD_TESTS_BASE<TV> BASE;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;

public:
    using BASE::last_frame;using BASE::grid;
    using BASE::particles;using BASE::frame_title;
    using BASE::write_substeps_level;using BASE::particles_per_cell;
    using BASE::substeps_delay_frame;using BASE::scale_mass;using BASE::data_directory;
    using BASE::output_directory;using BASE::restart;using BASE::dt;using BASE::stream_type;
    using BASE::time;using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;
    using BASE::order;using BASE::ghost;using BASE::use_affine;
    using BASE::cfl;using BASE::solver_tolerance;using BASE::Add_Source;
    using BASE::solver_iterations;using BASE::threads;using BASE::test_number;
    using BASE::resolution;using BASE::Seed_Particles;using BASE::Add_Particle;
    using BASE::regular_seeding;using BASE::Set_Grid;
    using BASE::random;using BASE::collision_objects;using BASE::user_resolution;
    using BASE::Add_Collision_Object;using BASE::Add_Fluid_Wall;using BASE::no_regular_seeding;
    using BASE::m;using BASE::s;using BASE::kg;using BASE::unit_p;using BASE::unit_mu;using BASE::unit_rho;
    using BASE::forced_collision_type;using BASE::destroy;using BASE::write_output_files;using BASE::read_output_files;
    using BASE::begin_frame;using BASE::end_frame;using BASE::density;
    using BASE::extra_T;using BASE::extra_int;using BASE::gravity;
    using BASE::side_bc_type;using BASE::BC_PERIODIC;using BASE::BC_FREE;using BASE::Test_dV;using BASE::bc_periodic;
    using BASE::use_analytic_field;using BASE::Seed_Particles_Analytic;
    using BASE::analytic_pressure;using BASE::Add_Pressure;
    using BASE::analytic_velocity;using BASE::Add_Velocity;using BASE::Check_Analytic_Velocity;
    using BASE::BC_SLIP;using BASE::BC_NOSLIP;using BASE::bc_velocity;using BASE::bc_pressure;
    using BASE::Setup_Analytic_Boundary_Conditions;
    using BASE::Add_Callbacks;using BASE::velocity;using BASE::mass;
    using BASE::viscosity;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS();

    void Initialize() override;

//#####################################################################
};
}

#endif
