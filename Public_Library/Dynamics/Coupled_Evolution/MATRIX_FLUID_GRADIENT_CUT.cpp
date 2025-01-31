//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_GRADIENT_CUT
//##################################################################### 
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Grids/SIDED_FACE_INDEX.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Dynamics/Coupled_Evolution/BOUNDARY_CONDITION_INFO.h>
#include <Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLECTION.h>
#include <Dynamics/Coupled_Evolution/MATRIX_FLUID_GRADIENT_CUT.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MATRIX_FLUID_GRADIENT_CUT<TV>::
MATRIX_FLUID_GRADIENT_CUT(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input)
    :BASE(index_map_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MATRIX_FLUID_GRADIENT_CUT<TV>::
~MATRIX_FLUID_GRADIENT_CUT()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MATRIX_FLUID_GRADIENT_CUT<TV>::
Compute(const ARRAY<bool,FACE_INDEX<d> >& psi_N_domain_boundary)
{
    // Will be set up by FLUID_TO_SOLID_INTERPOLATION_CUT or FLUID_TO_SOLID_INTERPOLATION_PHI
}
namespace PhysBAM{
template class MATRIX_FLUID_GRADIENT_CUT<VECTOR<float,1> >;
template class MATRIX_FLUID_GRADIENT_CUT<VECTOR<float,2> >;
template class MATRIX_FLUID_GRADIENT_CUT<VECTOR<float,3> >;
template class MATRIX_FLUID_GRADIENT_CUT<VECTOR<double,1> >;
template class MATRIX_FLUID_GRADIENT_CUT<VECTOR<double,2> >;
template class MATRIX_FLUID_GRADIENT_CUT<VECTOR<double,3> >;
}
