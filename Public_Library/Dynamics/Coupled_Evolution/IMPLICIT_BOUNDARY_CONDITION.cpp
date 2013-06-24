//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
using namespace PhysBAM;
template<class TV> IMPLICIT_BOUNDARY_CONDITION<TV>::
IMPLICIT_BOUNDARY_CONDITION()
{
}
template<class TV> IMPLICIT_BOUNDARY_CONDITION<TV>::
~IMPLICIT_BOUNDARY_CONDITION()
{
}
namespace PhysBAM{
template class IMPLICIT_BOUNDARY_CONDITION<VECTOR<float,1> >;
template class IMPLICIT_BOUNDARY_CONDITION<VECTOR<float,2> >;
template class IMPLICIT_BOUNDARY_CONDITION<VECTOR<float,3> >;
template class IMPLICIT_BOUNDARY_CONDITION<VECTOR<double,1> >;
template class IMPLICIT_BOUNDARY_CONDITION<VECTOR<double,2> >;
template class IMPLICIT_BOUNDARY_CONDITION<VECTOR<double,3> >;
}
