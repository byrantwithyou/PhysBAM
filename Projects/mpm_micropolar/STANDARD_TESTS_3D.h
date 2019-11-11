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
    using BASE::print_stats;using BASE::flip;using BASE::cfl;
    using BASE::mass;using BASE::gather_scatter;
    using BASE::test_diff;using BASE::threads;using BASE::test_number;
    using BASE::resolution;using BASE::Seed_Particles;using BASE::Add_Gravity;
    using BASE::Add_Quasi_Pressure;using BASE::random;
    using BASE::Add_Force;using BASE::data_directory;
    using BASE::stream_type;using BASE::force_helper;
    using BASE::Add_Particle;using BASE::quad_F_coeff;
    using BASE::m;using BASE::s;using BASE::kg;using BASE::unit_p;using BASE::unit_mu;using BASE::unit_rho;
    using BASE::scale_speed;using BASE::destroy;
    using BASE::write_output_files;using BASE::read_output_files;using BASE::begin_frame;
    using BASE::end_frame;using BASE::begin_time_step;using BASE::end_time_step;
    using BASE::Perturb;using BASE::Uniform;using BASE::extra_T;using BASE::extra_int;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS();

    void Write_Output_Files() override;
    void Read_Output_Files() override;
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
