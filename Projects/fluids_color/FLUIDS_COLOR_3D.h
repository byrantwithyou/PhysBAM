//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR_3D__
#define __FLUIDS_COLOR_3D__

#include "FLUIDS_COLOR_BASE.h"

namespace PhysBAM{

template<class T>
class FLUIDS_COLOR<VECTOR<T,3> >:public FLUIDS_COLOR_BASE<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    typedef FLUIDS_COLOR_BASE<TV> BASE;

public:
    using BASE::grid;using BASE::viewer_dir;using BASE::face_velocities;
    using BASE::write_substeps_level;using BASE::restart;using BASE::last_frame;using BASE::use_level_set_method;using BASE::use_pls;
    using BASE::dt;using BASE::levelset_color;using BASE::mu;using BASE::rho;using BASE::dump_matrix;using BASE::number_of_colors;
    using BASE::use_advection;using BASE::use_reduced_advection;using BASE::omit_solve;using BASE::use_discontinuous_velocity;
    using BASE::time_steps_per_frame;using BASE::use_p_null_mode;using BASE::Fill_Levelsets_From_Levelset_Color;
    using BASE::particle_levelset_evolution_multiple;using BASE::test_number;using BASE::use_pls_over_levelset;
    using BASE::analytic_velocity;using BASE::m;using BASE::s;using BASE::kg;using BASE::resolution;using BASE::analytic_levelset;
    using BASE::Large_Phi;using BASE::mu0;using BASE::mu1;using BASE::rho0;using BASE::rho1;using BASE::bc_type;using BASE::SLIP;
    using BASE::DIRICHLET;using BASE::NEUMANN;using BASE::surface_tension;using BASE::override_rho0;
    using BASE::override_rho1;using BASE::override_mu0;using BASE::override_mu1;using BASE::test_analytic_diff;
    using BASE::analytic_initial_only;using BASE::Set_Level_Set;using BASE::Level_Set_Error;
    using BASE::Velocity_Error;using BASE::Initialize_Common_Example;using BASE::After_Initialize_Example;

    FLUIDS_COLOR(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    ~FLUIDS_COLOR();

    void Initialize_Example();
//#####################################################################
};
}

#endif
