//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLAT_SYSTEM_FEM__
#define __FLAT_SYSTEM_FEM__
#include <Core/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV> class FLUID_LAYOUT_FEM;
template<class T,class TV>
void Compute_Full_Matrix(ARRAY<VECTOR<int,3> >& coded_entries,
    ARRAY<T>& code_values,ARRAY<T>& rhs_vector,FLUID_LAYOUT_FEM<TV>& fl,T mu,T unit_length);
}
#endif
