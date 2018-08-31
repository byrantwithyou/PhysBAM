//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLAT_SYSTEM_FEM__
#define __FLAT_SYSTEM_FEM__
#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV> struct FLUID_LAYOUT_FEM;
template<class TV> struct PARSE_DATA_FEM;
template<class TV> struct ANALYTIC_VECTOR;
template<class TV> struct ANALYTIC_SCALAR;

template<class T,class TV>
void Generate_Discretization(ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    ARRAY<T,CODE_ID>& code_values,const FLUID_LAYOUT_FEM<TV>& fl,
    const PARSE_DATA_FEM<TV>& pd,T mu,ARRAY<T,DOF_ID>& rhs,
    ANALYTIC_VECTOR<TV> * analytic_velocity,ANALYTIC_SCALAR<TV> * analytic_pressure);
}
#endif
