//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// #include <Tools/Grids_Uniform/CELL_ITERATOR.h>
// #include <Tools/Grids_Uniform/FACE_ITERATOR.h>
// #include <Tools/Log/DEBUG_SUBSTEPS.h>
// #include <Tools/Matrices/MATRIX.h>
// #include <Tools/Parsing/PARSE_ARGS.h>
// #include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
// #include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_CONST.h>
// #include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_LINE.h>
// #include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_NEST.h>
// #include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_ROTATE.h>
// #include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SCALE.h>
// #include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SIGNED.h>
// #include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_SPHERE.h>
// #include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_TRANSLATE.h>
// #include <Geometry/Analytic_Tests/ANALYTIC_LEVELSET_VORTEX.h>
// #include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
// #include <Geometry/Basic_Geometry/CYLINDER.h>
// #include <Geometry/Basic_Geometry/LINE_2D.h>
// #include <Geometry/Basic_Geometry/SPHERE.h>
// #include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
// #include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
// #include <Incompressible/Forces/VORTICITY_CONFINEMENT.h>
// #include <Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
// #include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include "FLUIDS_COLOR_3D.h"
// #ifdef USE_OPENMP
// #include <omp.h>
// #endif

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> FLUIDS_COLOR<VECTOR<T,3> >::
FLUIDS_COLOR(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :FLUIDS_COLOR_BASE<TV>(stream_type_input,parse_args)
{
    parse_args.Parse();

    if(!Initialize_Common_Example())
        Initialize_Example();

    After_Initialize_Example();
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> FLUIDS_COLOR<VECTOR<T,3> >::
~FLUIDS_COLOR()
{
}
//#####################################################################
// Function Initialize_Example
//#####################################################################
template<class T> void FLUIDS_COLOR<VECTOR<T,3> >::
Initialize_Example()
{
    switch(test_number){
        default: PHYSBAM_FATAL_ERROR("Missing test number");}
}
//#####################################################################
template class FLUIDS_COLOR<VECTOR<float,3> >;
template class FLUIDS_COLOR<VECTOR<double,3> >;
}
