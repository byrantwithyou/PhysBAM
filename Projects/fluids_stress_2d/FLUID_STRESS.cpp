//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
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
#include "FLUID_STRESS.h"

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> FLUID_STRESS<VECTOR<T,2> >::
FLUID_STRESS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :FLUID_STRESS_BASE<TV>(stream_type_input,parse_args),epsilon((T).1),radius((T).05),mode(2)
{
    parse_args.Parse();
    
    if(!Initialize_Common_Example())
        Initialize_Example();
    
    After_Initialize_Example();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FLUID_STRESS<VECTOR<T,2> >::
~FLUID_STRESS()
{
}
//#####################################################################
// Function Initialize_Example
//#####################################################################
template<class T> void FLUID_STRESS<VECTOR<T,2> >::
Initialize_Example()
{
    switch(test_number){
        case 100:
            grid.Initialize(TV_INT()+resolution,RANGE<TV>::Unit_Box()*(2*(T)pi)*m,true);
            analytic_levelset=new ANALYTIC_LEVELSET_CONST<TV>(-Large_Phi(),0,0);
            analytic_velocity=new ANALYTIC_VELOCITY_VORTEX<TV>();
            analytic_polymer_stress=new ANALYTIC_POLYMER_STRESS_CONST<TV>();
            break;
        default: PHYSBAM_FATAL_ERROR("Missing test number");}
}
//#####################################################################
template class FLUID_STRESS<VECTOR<float,2> >;
template class FLUID_STRESS<VECTOR<double,2> >;
}
