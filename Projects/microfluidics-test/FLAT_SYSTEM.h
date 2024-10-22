//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLAT_SYSTEM__
#define __FLAT_SYSTEM__
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR.h>
#include <functional>

namespace PhysBAM{

template<class T> class SYSTEM_MATRIX_HELPER;
template<class TV> class FLUID_LAYOUT;

template<class T,class TV>
void Compute_Full_Matrix(const GRID<TV>& grid,ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> >& coded_entries,
    ARRAY<T,CODE_ID>& code_values,ARRAY<T,DOF_ID>& rhs_vector,const FLUID_LAYOUT<TV>& fl,T mu);

template<class T,class TV>
void Solve_And_Display_Solution(const GRID<TV>& grid,const FLUID_LAYOUT<TV>& fl,
    const SYSTEM_MATRIX_HELPER<T>& MH,const ARRAY<T,DOF_ID>& rhs_vector,ARRAY<T,DOF_ID>* sol_vector,
    bool use_krylov,bool print_system,bool quiet);

}
#endif
