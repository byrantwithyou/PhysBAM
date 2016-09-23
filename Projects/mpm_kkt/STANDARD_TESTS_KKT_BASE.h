//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_KKT_BASE__
#define __STANDARD_TESTS_KKT_BASE__

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Deformables/Standard_Tests/DEFORMABLES_STANDARD_TESTS.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_KKT_EXAMPLE.h>
#include <functional>

namespace PhysBAM{

template<class TV> class STANDARD_TESTS_KKT;

template<class TV>
class STANDARD_TESTS_KKT_BASE:public MPM_KKT_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef MPM_KKT_EXAMPLE<TV> BASE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_VOLUME;

public:
    using BASE::initial_time;using BASE::last_frame;using BASE::grid;using BASE::particles;
    using BASE::mass;using BASE::force_helper;using BASE::use_oldroyd;
    using BASE::frame_title;using BASE::write_substeps_level;using BASE::gather_scatter;
    using BASE::collision_objects;using BASE::substeps_delay_frame;
    using BASE::output_directory;using BASE::mass_contour;
    using BASE::restart;using BASE::dt;using BASE::time;using BASE::use_early_gradient_transfer;
    using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;
    using BASE::ghost;using BASE::use_reduced_rasterization;using BASE::use_affine;;using BASE::use_f2p;
    using BASE::use_midpoint;using BASE::use_symplectic_euler;
    using BASE::print_stats;using BASE::flip;using BASE::cfl;using BASE::newton_tolerance;
    using BASE::newton_iterations;using BASE::solver_tolerance;using BASE::solver_iterations;
    using BASE::test_diff;using BASE::threads;using BASE::weights;
    using BASE::lagrangian_forces;using BASE::stream_type;
    using BASE::Add_Force;using BASE::Set_Weights;using BASE::deformable_body_collection;
    using BASE::Add_Collision_Object;using typename BASE::COLLISION_TYPE;using BASE::data_directory;
    using BASE::Add_Fluid_Wall;using BASE::quad_F_coeff;using BASE::coarse_grid;using BASE::use_FEM_mass;

    int test_number;
    int resolution;
    bool user_resolution;
    int stored_last_frame;
    bool user_last_frame;
    int order;
    int seed;
    int particles_per_cell;
    bool regular_seeding;
    T scale_mass;
    T scale_E;
    T scale_speed;
    T penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length;
    T penalty_damping_stiffness;
    RANDOM_NUMBERS<T> random;
    DEFORMABLES_STANDARD_TESTS<TV> tests;

    STANDARD_TESTS_KKT_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_KKT_BASE();

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

    void Add_Penalty_Collision_Object(IMPLICIT_OBJECT<TV>* io);
    template<class OBJECT> typename enable_if<!is_pointer<OBJECT>::value>::type
    Add_Penalty_Collision_Object(const OBJECT& object)
    {Add_Penalty_Collision_Object(new ANALYTIC_IMPLICIT_OBJECT<OBJECT>(object));}

    void Add_Particle(const TV& X,std::function<TV(const TV&)> V,std::function<MATRIX<T,TV::m>(const TV&)> dV,
        const T mass,const T volume);
    int Add_Gravity(TV g);
    int Add_Fixed_Corotated(T E,T nu,ARRAY<int>* affected_particles=0,bool no_mu=false);
    int Add_Neo_Hookean(T E,T nu,ARRAY<int>* affected_particles=0);
    void Add_Walls(int flags,COLLISION_TYPE type,T friction,T inset,bool penalty); // -x +x -y +y [ -z +z ], as bit flags
//#####################################################################
};
}

#endif
