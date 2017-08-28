//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_BASE__
#define __STANDARD_TESTS_BASE__

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_VECTOR.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Deformables/Standard_Tests/DEFORMABLES_STANDARD_TESTS.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_MAC_EXAMPLE.h>
#include <functional>

namespace PhysBAM{

template<class TV> class STANDARD_TESTS;
template<class T,int d> class ISOTROPIC_CONSTITUTIVE_MODEL;
template<class TV> class POISSON_DISK;
template<class TV> struct SOURCE_PATH;

template<class TV>
class STANDARD_TESTS_BASE:public MPM_MAC_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef MPM_MAC_EXAMPLE<TV> BASE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_VOLUME;

public:
//    typedef typename MPM_PARTICLE_SOURCE<TV>::PATH PATH;
    using BASE::initial_time;using BASE::last_frame;using BASE::grid;using BASE::particles;
    using BASE::frame_title;using BASE::write_substeps_level;using BASE::phases;
    using BASE::collision_objects;using BASE::substeps_delay_frame;using BASE::output_directory;
    using BASE::restart;using BASE::dt;using BASE::time;using BASE::lag_Dp;
    using BASE::frame_dt;using BASE::min_dt;using BASE::max_dt;using BASE::ghost;
    using BASE::use_affine;using BASE::only_write_particles;
    using BASE::cfl;using BASE::solver_tolerance;using BASE::solver_iterations;
    using BASE::threads;using BASE::weights;using BASE::stream_type;using BASE::Set_Weights;
    using BASE::Add_Collision_Object;using BASE::use_particle_volumes;using BASE::move_mass_inside;
    using BASE::test_system;using BASE::print_matrix;using BASE::move_mass_inside_nearest;
    using typename BASE::COLLISION_TYPE;using BASE::data_directory;using BASE::Add_Fluid_Wall;
    using BASE::test_output_prefix;using BASE::use_test_output;using BASE::flip;
    using BASE::begin_frame;using BASE::end_frame;using BASE::periodic_test_shift;using BASE::use_periodic_test_shift;
    using BASE::random;using BASE::use_viscosity;using typename BASE::PHASE;
    using BASE::use_constant_density;using BASE::bc_type;using BASE::BC_PERIODIC;
    using BASE::Add_Callbacks;using BASE::Print_Particle_Stats;using BASE::Print_Grid_Stats;
    using BASE::bc_velocity;using BASE::bc_pressure;using BASE::BC_WALL;
    using BASE::extrap_type;using BASE::use_extrap;

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
    bool override_output_directory;
    T m,s,kg;
    T unit_p,unit_rho,unit_mu;
    int forced_collision_type;
    ARRAY<std::function<void ()> > destroy;
    ARRAY<T> extra_T;
    ARRAY<int> extra_int;
    bool dump_collision_objects;
    bool test_diff;
    bool bc_periodic;
    bool use_analytic_field;
    ARRAY<ANALYTIC_VECTOR<TV>*,PHASE_ID> analytic_velocity;
    ARRAY<ANALYTIC_SCALAR<TV>*,PHASE_ID> analytic_pressure;
    T mu;
    bool analyze_u_modes=false;
    int dump_modes_freq=1;
    T max_ke=0;
    POISSON_DISK<TV>& poisson_disk;

    STANDARD_TESTS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_BASE();

    template<class F>
    void Add_Velocity(F f)
    {
        analytic_velocity.Append(Make_Analytic_Vector<TV>(f));
    }

    template<class F>
    void Add_Pressure(F f)
    {
        analytic_pressure.Append(Make_Analytic_Scalar<TV>(f));
    }

    void Setup_Analytic_Boundary_Conditions();
    
    virtual TV Compute_Analytic_Force(PHASE_ID p,const TV& X,T time) const;

    void Seed_Particles_Poisson(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell);

    template<class T_OBJECT> typename enable_if<!is_base_of<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>::type
    Seed_Particles_Poisson(const T_OBJECT& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles_Poisson(obj,V,dV,density,particles_per_cell);}

    void Seed_Particles_Uniform(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,const GRID<TV>& seed_grid);

    template<class T_OBJECT> typename enable_if<!is_base_of<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>::type
    Seed_Particles_Uniform(const T_OBJECT& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,const GRID<TV>& seed_grid)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles_Uniform(obj,V,dV,density,seed_grid);}

    void Seed_Particles(IMPLICIT_OBJECT<TV>& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell);

    template<class T_OBJECT> typename enable_if<!is_base_of<IMPLICIT_OBJECT<TV>,T_OBJECT>::value>::type
    Seed_Particles(const T_OBJECT& object,std::function<TV(const TV&)> V,
        std::function<MATRIX<T,TV::m>(const TV&)> dV,T density,T particles_per_cell)
    {ANALYTIC_IMPLICIT_OBJECT<T_OBJECT> obj(object);Seed_Particles(obj,V,dV,density,particles_per_cell);}

    template<class T_OBJECT> auto
    Seed_Particles_Analytic(const T_OBJECT& object,PHASE_ID pid,T density,T particles_per_cell)
    {return Seed_Particles(object,
            [this,pid](const TV& X){return analytic_velocity(pid)->v(X,0);},
            [this,pid](const TV& X){return analytic_velocity(pid)->dX(X,0);},
            density,particles_per_cell);}

    void Add_Particle(const TV& X,std::function<TV(const TV&)> V,std::function<MATRIX<T,TV::m>(const TV&)> dV,
        const T mass,const T volume);
    void Add_Particle(const TV& X,const TV& V,const MATRIX<T,TV::m>& dV,const T mass,const T volume);

    void Set_Grid(const RANGE<TV>& domain,TV_INT resolution_scale=TV_INT()+1,int default_resolution=32);
    void Set_Grid(const RANGE<TV>& domain,TV_INT resolution_scale,TV_INT resolution_padding,
        int resolution_multiple=1,int default_resolution=32);
    void Test_dV(std::function<TV(const TV&)> V,std::function<MATRIX<T,TV::m>(const TV&)> dV) const;
    void Set_Phases(const ARRAY<T,PHASE_ID>& phase_densities);
    void Check_Analytic_Velocity() const;
    void Velocity_Fourier_Analysis() const;
    void Add_Source(const TV& X0,const TV& n,IMPLICIT_OBJECT<TV>* io,
        std::function<void(TV X,T ts,T t,SOURCE_PATH<TV>& p)> path,T density,
        T particles_per_cell,bool owns_io);
//#####################################################################
};
}

#endif
