//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_BASE__
#define __STANDARD_TESTS_BASE__

#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <boost/function.hpp>

namespace PhysBAM{

template<class TV> class STANDARD_TESTS;

template<class TV>
class STANDARD_TESTS_BASE:public MPM_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef MPM_EXAMPLE<TV> BASE;

public:
    using BASE::initial_time;using BASE::last_frame;using BASE::grid;using BASE::particles;
    using BASE::frame_title;using BASE::write_substeps_level;using BASE::gather_scatter;
    using BASE::substeps_delay_frame;using BASE::write_output_files;
    using BASE::output_directory;using BASE::restart;using BASE::dt;using BASE::time;
    using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;
    using BASE::ghost;using BASE::use_reduced_rasterization;using BASE::use_affine;
    using BASE::use_midpoint;using BASE::flip;using BASE::cfl;using BASE::newton_tolerance;
    using BASE::newton_iterations;using BASE::solver_tolerance;using BASE::solver_iterations;
    using BASE::test_diff;using BASE::threads;using BASE::weights;
    using BASE::Add_Force;

    int test_number;
    int resolution;
    int stored_last_frame;
    bool user_last_frame;
    int order;
    int seed;
    int particles_per_cell;
    RANDOM_NUMBERS<T> random;

    STANDARD_TESTS_BASE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_BASE();

    void Seed_Particles(IMPLICIT_OBJECT<TV>& object,boost::function<TV(const TV&)> V,
        boost::function<MATRIX<T,TV::m>(const TV&)> dV,T density,int particles_per_cell);

    template<class T_OBJECT> typename DISABLE_IF<IS_BASE_OF<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>::TYPE
    Seed_Particles(const T_OBJECT& object,boost::function<TV(const TV&)> V,
        boost::function<MATRIX<T,TV::m>(const TV&)> dV,T density,int particles_per_cell)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles(obj,V,dV,density,particles_per_cell);}

    int Add_Gravity(TV g);
    int Add_Fixed_Corotated(T E,T nu,ARRAY<int>* affected_particles=0);

//#####################################################################
};
}

#endif
