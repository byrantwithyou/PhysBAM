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
    using BASE::lagrangian_forces;using BASE::mass;using BASE::gather_scatter;
    using BASE::solver_iterations;using BASE::test_diff;using BASE::threads;using BASE::test_number;
    using BASE::resolution;using BASE::Seed_Particles;using BASE::Add_Gravity;
    using BASE::Add_Fixed_Corotated;using BASE::random;using BASE::Seed_Lagrangian_Particles;
    using BASE::Add_Force;using BASE::Add_Walls;using BASE::data_directory;
    using BASE::stream_type;using BASE::use_oldroyd;using BASE::force_helper;
    using BASE::Seed_Particles_Helper;using BASE::Add_Neo_Hookean;using BASE::Add_Collision_Object;
    using BASE::Add_Particle;using BASE::Add_Penalty_Collision_Object;using BASE::scale_speed;

    int Nsurface;
    ARRAY<int> steal;

    STANDARD_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS();

    void Write_Output_Files(const int frame) override;
    void Read_Output_Files(const int frame) override;
    void Initialize() override;
    void Begin_Frame(const int frame) override;
    void End_Frame(const int frame) override;
    void Begin_Time_Step(const T time) override;
    void End_Time_Step(const T time) override;

    LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* Initialize_Implicit_Surface(TRIANGULATED_SURFACE<T>& surface,int max_res);

    // additional storage
    int foo_int1;
    T foo_T1;
    T foo_T2;
    T foo_T3;
    TRIANGULATED_SURFACE<T>* foo_surface1;
    TRIANGULATED_SURFACE<T>* foo_surface2;
    LEVELSET_IMPLICIT_OBJECT<TV>* foo_levelset1;

//#####################################################################
};
}

#endif
