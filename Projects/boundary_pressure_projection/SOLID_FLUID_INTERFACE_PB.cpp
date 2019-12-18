//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include "SOLID_FLUID_INTERFACE_PB.h"
namespace PhysBAM{

//#####################################################################
// Function Compute_BC
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Compute_BC(FLUID_SOLVER<TV>* fluid_solver,SOLID_BC<TV>* solid_bc,T time,T dt) const
{
    LOG::printf("DO THIS: %s\n",__PRETTY_FUNCTION__);
}

//#####################################################################
// Function Compute_BC
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Compute_BC(SOLID_SOLVER<TV>* solid_solver,FLUID_BC<TV>* fluid_bc,T time,T dt) const
{
    LOG::printf("DO THIS: %s\n",__PRETTY_FUNCTION__);
}

//#####################################################################
// Function Interpolate_Velocity
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Interpolate_Velocity(FLUID_BOUNDARY_VECTOR<TV>* u, const SOLID_BOUNDARY_VECTOR<TV>* v)
{
}

//#####################################################################
// Function Distribute_Force
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Distribute_Force(SOLID_BOUNDARY_VECTOR<TV>* v, const FLUID_BOUNDARY_VECTOR<TV>* u)
{
}

//#####################################################################
// Function Get_Boundary
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Get_Boundary(SOLID_BOUNDARY_VECTOR<TV>* v)
{
}

template class SOLID_FLUID_INTERFACE_PB<VECTOR<float,2> >;
template class SOLID_FLUID_INTERFACE_PB<VECTOR<float,3> >;
template class SOLID_FLUID_INTERFACE_PB<VECTOR<double,2> >;
template class SOLID_FLUID_INTERFACE_PB<VECTOR<double,3> >;
}
