//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_BASE__
#define __STANDARD_TESTS_BASE__

#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include "PBD_EXAMPLE.h"

namespace PhysBAM{

template<class TV> class STANDARD_TESTS;

template<class TV>
class STANDARD_TESTS_BASE:public PBD_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef PBD_EXAMPLE<TV> BASE;

public:
    using BASE::initial_time;using BASE::last_frame;
    using BASE::frame_title;using BASE::write_substeps_level;
    using BASE::substeps_delay_frame;
    using BASE::output_directory;using BASE::restart;using BASE::dt;using BASE::time;
    using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;
    using BASE::print_stats;
    using BASE::solver_iterations;
    using BASE::test_diff;using BASE::threads;
    using BASE::X;
    using BASE::V;
    using BASE::w;

    int test_number;
    int resolution;
    bool user_resolution;
    int stored_last_frame;
    bool user_last_frame;
    int seed;
    T scale_mass;
    T scale_stiffness;
    T scale_speed;
    RANDOM_NUMBERS<T> random;

    STANDARD_TESTS_BASE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_BASE();
//#####################################################################
};
}

#endif
