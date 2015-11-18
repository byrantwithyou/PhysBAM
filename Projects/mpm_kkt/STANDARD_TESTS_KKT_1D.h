//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_KKT_1D__
#define __STANDARD_TESTS_KKT_1D__

#include <Tools/Parsing/PARSE_ARGS.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_KKT_EXAMPLE.h>
#include "STANDARD_TESTS_KKT_BASE.h"

namespace PhysBAM{

template<class TV> class STANDARD_TESTS_KKT;

template<class T>
class STANDARD_TESTS_KKT<VECTOR<T,1> >:public STANDARD_TESTS_KKT_BASE<VECTOR<T,1> >
{
    typedef VECTOR<T,1> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef STANDARD_TESTS_KKT_BASE<TV> BASE;
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;

public:
    using BASE::initial_time;using BASE::last_frame;using BASE::grid;using BASE::particles;
    using BASE::mass;using BASE::force_helper;
    using BASE::frame_title;using BASE::write_substeps_level;using BASE::particles_per_cell;
    using BASE::substeps_delay_frame;using BASE::scale_mass;using BASE::scale_E;using BASE::scale_speed;
    using BASE::output_directory;using BASE::restart;using BASE::dt;using BASE::time;
    using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;using BASE::order;
    using BASE::ghost;using BASE::use_reduced_rasterization;using BASE::use_affine;
    using BASE::use_midpoint;using BASE::print_stats;using BASE::flip;using BASE::cfl;
    using BASE::newton_tolerance;using BASE::newton_iterations;using BASE::solver_tolerance;
    using BASE::solver_iterations;using BASE::test_diff;using BASE::threads;
    using BASE::lagrangian_forces;using BASE::use_oldroyd;using BASE::penalty_collisions_stiffness;
    using BASE::test_number;using BASE::resolution;using BASE::Seed_Particles;using BASE::Add_Gravity;
    using BASE::Add_Particle;using BASE::gather_scatter;using BASE::penalty_damping_stiffness;
    using BASE::Add_Fixed_Corotated;using BASE::Add_Neo_Hookean;using BASE::random;
    using BASE::collision_objects;using BASE::user_resolution;using BASE::Add_Walls;
    using BASE::tests;using BASE::Seed_Lagrangian_Particles;using BASE::Add_Collision_Object;
    using BASE::Add_Force;using BASE::Add_Fluid_Wall;using BASE::quad_F_coeff;using BASE::coarse_grid;

    // surface tension stuff
    bool use_surface_tension;
    int Nsurface;
    ARRAY<int> steal;

    STANDARD_TESTS_KKT(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_KKT();

    void Write_Output_Files(const int frame) override;
    void Read_Output_Files(const int frame) override;
    void Initialize() override;
    void Begin_Frame(const int frame) override;
    void End_Frame(const int frame) override;
    void Begin_Time_Step(const T time) override;
    void End_Time_Step(const T time) override;

//#####################################################################
};
}

#endif
