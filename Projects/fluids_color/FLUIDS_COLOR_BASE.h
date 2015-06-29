//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUIDS_COLOR_BASE__
#define __FLUIDS_COLOR_BASE__

#include <Tools/Parsing/PARSE_ARGS.h>
#include <Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
#include "ANALYTIC_POLYMER_STRESS.h"
#include "ANALYTIC_VELOCITY.h"

namespace PhysBAM{

template<class TV> class FLUIDS_COLOR;

template<class TV>
class FLUIDS_COLOR_BASE:public PLS_FC_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef PLS_FC_EXAMPLE<TV> BASE;
    typedef SYMMETRIC_MATRIX<T,TV::m> SM;

public:
    using BASE::grid;using BASE::output_directory;using BASE::face_velocities;using BASE::write_substeps_level;
    using BASE::restart;using BASE::last_frame;using BASE::use_level_set_method;using BASE::use_pls;
    using BASE::dt;using BASE::levelset_color;using BASE::mu;using BASE::rho;using BASE::dump_matrix;
    using BASE::sparse_dump_matrix;using BASE::number_of_colors;using BASE::num_multigrid_levels;using BASE::use_multigrid;
    using BASE::use_advection;using BASE::use_reduced_advection;using BASE::omit_solve;using BASE::use_discontinuous_velocity;
    using BASE::time_steps_per_frame;using BASE::use_p_null_mode;using BASE::Fill_Levelsets_From_Levelset_Color;
    using BASE::particle_levelset_evolution_multiple;using BASE::face_color;using BASE::substeps_delay_frame;
    using BASE::dump_largest_eigenvector;using BASE::save_pressure;using BASE::use_polymer_stress;using BASE::pressure;
    using BASE::polymer_stress;using BASE::polymer_stress_coefficient;using BASE::inv_Wi;using BASE::max_iter;
    using BASE::test_system;using BASE::use_preconditioner;using BASE::solver_tolerance;

    enum WORKAROUND{SLIP=-3,DIRICHLET=-2,NEUMANN=-1}; // From CELL_DOMAIN_INTERFACE_COLOR

    int test_number;
    int resolution;
    int stored_last_frame;
    bool user_last_frame;
    T mu0,mu1;
    T rho0,rho1;
    T inv_Wi0,inv_Wi1;
    T beta0,beta1;
    T unit_mu,unit_rho,unit_st,unit_p;
    T weiss,weiss_inv;
    T m,s,kg;
    int bc_type;
    bool bc_n,bc_d,bc_s;
    bool test_analytic_diff;
    int refine;
    static T Large_Phi() {return 1000;}
    T surface_tension;
    bool override_rho0;
    bool override_rho1;
    bool override_mu0;
    bool override_mu1;
    bool override_inv_Wi0,override_inv_Wi1;
    bool override_beta0,override_beta1;
    bool override_surface_tension;
    bool use_pls_over_levelset;
    bool use_levelset_over_pls;
    
    TV gravity;

    ARRAY<ANALYTIC_VELOCITY<TV>*> analytic_velocity,initial_analytic_velocity;
    ARRAY<ANALYTIC_PRESSURE<TV>*> analytic_pressure,initial_analytic_pressure;
    ARRAY<ANALYTIC_POLYMER_STRESS<TV>*> analytic_polymer_stress;
    ANALYTIC_LEVELSET<TV>* analytic_levelset;
    bool analytic_initial_only;
    int number_of_threads;
    bool override_output_directory;

    bool use_u0,use_u1,use_p0,use_p1,use_S0,use_S1;
    std::string str_u0,str_u1,str_p0,str_p1,str_S0,str_S1;

    FLUIDS_COLOR_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    ~FLUIDS_COLOR_BASE();

    template<class F>
    void Add_Velocity(F f)
    {
        analytic_velocity.Append(Make_Velocity<TV>(f));
    }

    template<class F>
    void Add_Pressure(F f)
    {
        analytic_pressure.Append(Make_Pressure<TV>(f));
    }

    template<class F>
    void Add_Polymer_Stress(F f)
    {
        analytic_polymer_stress.Append(Make_Polymer_Stress<TV>(f));
    }

    void After_Initialize_Example();
    bool Initialize_Common_Example();
    void Write_Output_Files(const int frame);
    void Initialize();
    void Set_Level_Set(T time);
    void Level_Set_Error(T time);
    void Velocity_Error(T time);
    void Dump_Analytic_Levelset(T time);
    void Get_Initial_Velocities(T time);
    void Get_Initial_Polymer_Stresses(T time);
    void Analytic_Test();
    void Begin_Time_Step(const T time) override;
    void End_Time_Step(const T time) override;
    MATRIX<T,TV::m> Stress(const TV& X,int color,T time);
    SYMMETRIC_MATRIX<T,TV::m> Polymer_Stress(const TV& X,int color,T time);
    SYMMETRIC_MATRIX<T,TV::m> Polymer_Stress_Forcing_Term(const TV& X,int color,T time);
    TV Jump_Interface_Condition(const TV& X,int color0,int color1,T time) override;
    TV Volume_Force(const TV& X,int color,T time) override;
    TV Velocity_Jump(const TV& X,int color0,int color1,T time) override;
//#####################################################################
};
}
#endif
