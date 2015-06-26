//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR_2D__
#define __FLUIDS_COLOR_2D__

#include "FLUIDS_COLOR_BASE.h"

namespace PhysBAM{

template<class T>
class FLUIDS_COLOR<VECTOR<T,2> >:public FLUIDS_COLOR_BASE<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    typedef FLUIDS_COLOR_BASE<TV> BASE;
    typedef SYMMETRIC_MATRIX<T,TV::m> SM;

public:
    using BASE::grid;using BASE::use_level_set_method;using BASE::use_p_null_mode;using BASE::test_number;
    using BASE::analytic_velocity;using BASE::m;using BASE::s;using BASE::kg;using BASE::resolution;
    using BASE::analytic_levelset;using BASE::Large_Phi;using BASE::mu0;using BASE::mu1;using BASE::rho0;
    using BASE::rho1;using BASE::bc_type;using BASE::SLIP;using BASE::DIRICHLET;using BASE::NEUMANN;
    using BASE::unit_rho;using BASE::unit_mu;using BASE::unit_st;using BASE::surface_tension;
    using BASE::override_rho0;using BASE::override_rho1;using BASE::override_mu0;using BASE::override_mu1;
    using BASE::override_inv_Wi0;using BASE::override_beta0;
    using BASE::test_analytic_diff;using BASE::Initialize_Common_Example;using BASE::After_Initialize_Example;
    using BASE::use_discontinuous_velocity;using BASE::gravity;using BASE::analytic_initial_only;
    using BASE::override_surface_tension;using BASE::unit_p;using BASE::use_advection;
    using BASE::use_polymer_stress;using BASE::analytic_polymer_stress;using BASE::Add_Polymer_Stress;
    using BASE::polymer_stress_coefficient;using BASE::inv_Wi;using BASE::Add_Velocity;using BASE::Add_Pressure;

    T epsilon,radius;
    int mode;

    FLUIDS_COLOR(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    ~FLUIDS_COLOR();

    void Initialize_Example();
    void Add_Vortex(T mu,T rho,TV trans_vel=TV(),T new_l=0);
//#####################################################################
};
}

#endif
