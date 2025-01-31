//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Incompressible/Solids_And_Fluids/BOUNDARY_CONDITIONS_CALLBACKS.h>
#include <Dynamics/Coupled_Evolution/EXAMPLE_BOUNDARY_CONDITION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> EXAMPLE_BOUNDARY_CONDITION<TV>::
EXAMPLE_BOUNDARY_CONDITION(BOUNDARY_CONDITIONS_CALLBACKS<TV>* callback_input)
    :callback(callback_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> EXAMPLE_BOUNDARY_CONDITION<TV>::
~EXAMPLE_BOUNDARY_CONDITION()
{
}
//#####################################################################
// Function Update_Boundary_Conditions
//#####################################################################
template<class TV> void EXAMPLE_BOUNDARY_CONDITION<TV>::
Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& p,
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time)
{
    callback->Mark_Outside(psi_D);
}
namespace PhysBAM{
template class EXAMPLE_BOUNDARY_CONDITION<VECTOR<float,1> >;
template class EXAMPLE_BOUNDARY_CONDITION<VECTOR<float,2> >;
template class EXAMPLE_BOUNDARY_CONDITION<VECTOR<float,3> >;
template class EXAMPLE_BOUNDARY_CONDITION<VECTOR<double,1> >;
template class EXAMPLE_BOUNDARY_CONDITION<VECTOR<double,2> >;
template class EXAMPLE_BOUNDARY_CONDITION<VECTOR<double,3> >;
}
