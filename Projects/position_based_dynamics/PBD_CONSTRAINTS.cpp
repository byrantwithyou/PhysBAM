//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "PBD_CONSTRAINTS.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PBD_CONSTRAINTS_BASE<TV>::
PBD_CONSTRAINTS_BASE()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PBD_CONSTRAINTS_BASE<TV>::
~PBD_CONSTRAINTS_BASE()
{
}
template class PBD_CONSTRAINTS_BASE<VECTOR<float,2> >;
template class PBD_CONSTRAINTS_BASE<VECTOR<double,2> >;
template class PBD_CONSTRAINTS_BASE<VECTOR<float,3> >;
template class PBD_CONSTRAINTS_BASE<VECTOR<double,3> >;
}
