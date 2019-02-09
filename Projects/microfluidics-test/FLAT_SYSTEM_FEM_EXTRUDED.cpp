//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include "COMMON.h"
#include "FLUID_LAYOUT_FEM_EXTRUDED.h"
#include "FLAT_SYSTEM_FEM_EXTRUDED.h"

namespace PhysBAM{

template<class T,class TV,class TV2>
void Generate_Discretization(ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    ARRAY<T,CODE_ID>& code_values,const FLUID_LAYOUT_FEM<TV>& fl,
    const PARSE_DATA_FEM<TV2,TV>& pd,T mu,ARRAY<T,DOF_ID>& rhs)
{
}

template<class T,class TV,class TV2>
void Solve_And_Display_Solution(const FLUID_LAYOUT_FEM<TV>& fl,const PARSE_DATA_FEM<TV2,TV>& pd,
    const SYSTEM_MATRIX_HELPER<T>& MH,const ARRAY<T,DOF_ID>& rhs_vector,
    ARRAY<T,DOF_ID>* sol_out)
{
}

template void Generate_Discretization<double,VECTOR<double,3>,VECTOR<double,2> >(
    ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID>,int>&,ARRAY<double,CODE_ID>&,FLUID_LAYOUT_FEM<VECTOR<double,3> > const&,
    PARSE_DATA_FEM<VECTOR<double,2>,VECTOR<double,3> > const&,double,ARRAY<double,DOF_ID>&);
template void Solve_And_Display_Solution<double,VECTOR<double,3>,VECTOR<double,2> >(
    FLUID_LAYOUT_FEM<VECTOR<double,3> > const&,PARSE_DATA_FEM<VECTOR<double,2>,VECTOR<double,3> > const&,
    SYSTEM_MATRIX_HELPER<double> const&,ARRAY<double,DOF_ID> const&,ARRAY<double,DOF_ID>*);
}
