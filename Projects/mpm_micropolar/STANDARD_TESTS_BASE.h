//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_BASE__
#define __STANDARD_TESTS_BASE__

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MICROPOLAR_EXAMPLE.h>
#include <functional>

namespace PhysBAM{

template<class TV> class STANDARD_TESTS;
template<class T,int d> class ISOTROPIC_CONSTITUTIVE_MODEL;

template<class TV>
class STANDARD_TESTS_BASE:public MPM_MICROPOLAR_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef MPM_MICROPOLAR_EXAMPLE<TV> BASE;

public:
    using BASE::last_frame;using BASE::grid;using BASE::particles;
    using BASE::mass;using BASE::force_helper;
    using BASE::frame_title;using BASE::write_substeps_level;using BASE::gather_scatter;
    using BASE::substeps_delay_frame;
    using BASE::viewer_dir;
    using BASE::restart;using BASE::dt;using BASE::time;
    using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;
    using BASE::ghost;using BASE::use_affine;using BASE::cfl_F;using BASE::use_strong_cfl;
    using BASE::test_output_prefix;
    using BASE::print_stats;using BASE::only_write_particles;using BASE::use_test_output;
    using BASE::flip;using BASE::cfl;
    using BASE::test_diff;using BASE::threads;using BASE::weights;
    using BASE::stream_type;
    using BASE::Add_Force;using BASE::Set_Weights;
    using BASE::data_directory;using BASE::reflection_bc_flags;
    using BASE::quad_F_coeff;using BASE::use_sound_speed_cfl;using BASE::cfl_sound;

    int test_number;
    int resolution;
    bool user_resolution;
    int stored_last_frame;
    bool user_last_frame;
    int order;
    int seed;
    T particles_per_cell;
    bool regular_seeding,no_regular_seeding;
    T scale_mass;
    T scale_E;
    T scale_speed;
    bool override_output_directory;
    T m,s,kg;
    T unit_p,unit_rho,unit_mu;
    std::function<void ()> write_output_files;
    std::function<void ()> read_output_files;
    std::function<void ()> destroy;
    ARRAY<T> extra_T;
    ARRAY<int> extra_int;

    RANDOM_NUMBERS<T> random;

    STANDARD_TESTS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_BASE();

    T Perturb(T a);
    T Uniform(T a,T b);

    void Seed_Particles_Poisson(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell);

    template<class T_OBJECT> enable_if_t<!is_base_of<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>
    Seed_Particles_Poisson(const T_OBJECT& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles_Poisson(obj,V,dV,density,particles_per_cell);}

    void Seed_Particles_Uniform(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,const GRID<TV>& seed_grid);

    template<class T_OBJECT> enable_if_t<!is_base_of<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>
    Seed_Particles_Uniform(const T_OBJECT& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,const GRID<TV>& seed_grid)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles_Uniform(obj,V,dV,density,seed_grid);}

    void Seed_Particles(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell);

    template<class T_OBJECT> enable_if_t<!is_base_of<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>
    Seed_Particles(const T_OBJECT& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles(obj,V,dV,density,particles_per_cell);}

    void Add_Particle(const TV& X,std::function<TV(const TV&)> V,std::function<MATRIX<T,TV::m>(const TV&)> dV,
        const T mass,const T volume);
    int Add_Gravity(TV g);
    int Add_Quasi_Pressure(T E,T gamma);
    void Set_Grid(const RANGE<TV>& domain,TV_INT resolution_scale=TV_INT()+1,int default_resolution=32);
    void Set_Grid(const RANGE<TV>& domain,TV_INT resolution_scale,TV_INT resolution_padding,
        int resolution_multiple=1,int default_resolution=32);
//#####################################################################
};
}

#endif
