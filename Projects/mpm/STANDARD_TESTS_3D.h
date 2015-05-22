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

template<class T>
class STANDARD_TESTS<VECTOR<T,3> >:public STANDARD_TESTS_BASE<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef STANDARD_TESTS_BASE<TV> BASE;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;

public:
    using BASE::initial_time;using BASE::last_frame;using BASE::grid;using BASE::particles;
    using BASE::frame_title;using BASE::write_substeps_level;using BASE::particles_per_cell;
    using BASE::substeps_delay_frame;using BASE::scale_mass;using BASE::scale_E;
    using BASE::output_directory;using BASE::restart;using BASE::dt;using BASE::time;
    using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;using BASE::order;
    using BASE::ghost;using BASE::use_reduced_rasterization;using BASE::use_affine;
    using BASE::use_midpoint;using BASE::print_stats;using BASE::flip;using BASE::cfl;
    using BASE::newton_tolerance;using BASE::newton_iterations;using BASE::solver_tolerance;
    using BASE::solver_iterations;using BASE::test_diff;using BASE::threads;using BASE::test_number;
    using BASE::resolution;using BASE::Seed_Particles;using BASE::Add_Gravity;
    using BASE::Add_Fixed_Corotated;using BASE::random;using BASE::Seed_Lagrangian_Particles;
    using BASE::Add_Force;using BASE::Add_Walls;using BASE::data_directory;

    STANDARD_TESTS(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS();

    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
    void Read_Output_Files(const int frame) PHYSBAM_OVERRIDE;
    void Initialize() PHYSBAM_OVERRIDE;
    void Begin_Frame(const int frame) PHYSBAM_OVERRIDE;
    void End_Frame(const int frame) PHYSBAM_OVERRIDE;
    void Begin_Time_Step(const T time) PHYSBAM_OVERRIDE;
    void End_Time_Step(const T time) PHYSBAM_OVERRIDE;

//#####################################################################
};
}

#endif
