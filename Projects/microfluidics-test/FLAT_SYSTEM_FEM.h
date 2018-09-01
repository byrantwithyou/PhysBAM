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

template<class T,class TV>
void Generate_Discretization(ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    ARRAY<T,CODE_ID>& code_values,const FLUID_LAYOUT_FEM<TV>& fl,
    const PARSE_DATA_FEM<TV>& pd,T mu,ARRAY<T,DOF_ID>& rhs);
template<class T,class TV>
void Solve_And_Display_Solution(const FLUID_LAYOUT_FEM<TV>& fl,const PARSE_DATA_FEM<TV>& pd,
    const SYSTEM_MATRIX_HELPER<T>& MH,const ARRAY<T,DOF_ID>& rhs_vector,
    ARRAY<T,DOF_ID>* sol_out=0);
}
#endif
