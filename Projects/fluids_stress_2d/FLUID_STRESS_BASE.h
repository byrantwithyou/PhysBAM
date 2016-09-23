//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_STRESS_BASE__
#define __FLUID_STRESS_BASE__

#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_BOX.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ELLIPSOID.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_LINE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_NEST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ROTATE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SCALE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SPHERE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_TRANSLATE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_VORTEX.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include "ANALYTIC_POLYMER_STRESS.h"
#include "ANALYTIC_VELOCITY.h"
#include "STRESS_EXAMPLE.h"

namespace PhysBAM{

template<class TV> class FLUID_STRESS;

template<class TV>
class FLUID_STRESS_BASE:public STRESS_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef STRESS_EXAMPLE<TV> BASE;

public:
    using BASE::grid;using BASE::output_directory;using BASE::face_velocities;using BASE::write_substeps_level;
    using BASE::restart;using BASE::last_frame;using BASE::dt;using BASE::levelset;
    using BASE::time_steps_per_frame;using BASE::substeps_delay_frame;using BASE::polymer_stress;
    using BASE::number_of_ghost_cells;using BASE::inv_Wi;using BASE::use_du_terms;

    enum WORKAROUND{SLIP=-3,DIRICHLET=-2,NEUMANN=-1}; // From CELL_DOMAIN_INTERFACE_COLOR

    int test_number;
    int resolution;
    int parameter;
    int stored_last_frame;
    bool user_last_frame;
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
    bool override_surface_tension;
    bool use_pls_over_levelset;
    bool use_levelset_over_pls;
    bool use_inv_Wi;

    TV gravity;

    ANALYTIC_VELOCITY<TV>* analytic_velocity;
    ANALYTIC_POLYMER_STRESS<TV>* analytic_polymer_stress;
    ANALYTIC_LEVELSET_SIGNED<TV>* analytic_levelset;
    bool analytic_initial_only;
    int number_of_threads;
    bool override_output_directory;

    FLUID_STRESS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    ~FLUID_STRESS_BASE();

    void After_Initialize_Example();
    bool Initialize_Common_Example();
    void Write_Output_Files(const int frame) override;
    void Initialize() override;
    void Set_Level_Set(T time);
    void Stress_Error(T time);
    void Dump_Analytic_Levelset(T time);
    void Get_Velocities(T time) override;
    void Get_Initial_Polymer_Stresses() override;
    void Analytic_Test();
    void Begin_Time_Step(const T time) override;
    void End_Time_Step(const T time) override;
    MATRIX<T,TV::m> Stress(const TV& X,T time);
    SYMMETRIC_MATRIX<T,TV::m> Polymer_Stress(const TV& X,T time) override;
    SYMMETRIC_MATRIX<T,TV::m> Polymer_Stress_Forcing_Term(const TV& X,T time) override;
//#####################################################################
};
}

#endif
