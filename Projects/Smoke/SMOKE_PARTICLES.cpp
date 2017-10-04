//#####################################################################
// Copyright 2015, Peter Chen.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "SMOKE_PARTICLES.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SMOKE_PARTICLES<TV>::
SMOKE_PARTICLES()
{
    this->Store_Velocity();
    this->Store_Mass();
    Add_Array("smoke_C",&C);
    Add_Array("smoke_X0",&X0);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SMOKE_PARTICLES<TV>::
~SMOKE_PARTICLES()
{}
//#####################################################################
template class SMOKE_PARTICLES<VECTOR<float,2> >;
template class SMOKE_PARTICLES<VECTOR<float,3> >;
template class SMOKE_PARTICLES<VECTOR<double,2> >;
template class SMOKE_PARTICLES<VECTOR<double,3> >;
}
