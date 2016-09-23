//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLUID_STRESS__
#define __FLUID_STRESS__

#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_BOX.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ELLIPSOID.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_LINE.h>
#include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_MOVING_ELLIPSOID.h>
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
#include "FLUID_STRESS_BASE.h"

namespace PhysBAM{

template<class T>
class FLUID_STRESS<VECTOR<T,2> >:public FLUID_STRESS_BASE<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    typedef FLUID_STRESS_BASE<TV> BASE;

public:
    using BASE::grid;using BASE::test_number;
    using BASE::analytic_velocity;using BASE::m;using BASE::s;using BASE::kg;using BASE::resolution;
    using BASE::analytic_levelset;using BASE::Large_Phi;using BASE::bc_type;using BASE::SLIP;using BASE::DIRICHLET;using BASE::NEUMANN;
    using BASE::unit_rho;using BASE::unit_mu;using BASE::unit_st;using BASE::surface_tension;
    using BASE::override_rho0;using BASE::override_rho1;using BASE::override_mu0;using BASE::override_mu1;
    using BASE::test_analytic_diff;using BASE::Initialize_Common_Example;using BASE::After_Initialize_Example;
    using BASE::analytic_initial_only;using BASE::unit_p;using BASE::use_advection;
    using BASE::analytic_polymer_stress;

    T epsilon,radius;
    int mode;

    FLUID_STRESS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    ~FLUID_STRESS();

    void Initialize_Example();
//#####################################################################
};
}
#endif
