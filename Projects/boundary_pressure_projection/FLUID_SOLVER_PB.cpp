//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "FLUID_BC_PB.h"
#include "FLUID_SOLVER_PB.h"
#include "FLUID_STATE_PB.h"
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_SOLVER_PB<TV>::
FLUID_SOLVER_PB()
{
}

//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_SOLVER_PB<TV>::
~FLUID_SOLVER_PB()
{
}

//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Initialize()
{
    driver->Initialize();
}

//#####################################################################
// Function Write
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Write(int frame) const
{
    driver->Write_Output_Files(frame);
}

//#####################################################################
// Function Read
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Read(int frame)
{
}

//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> auto FLUID_SOLVER_PB<TV>::
Compute_Dt(T time) const -> T
{
    return driver->Compute_Fluids_Dt(time);
}

//#####################################################################
// Function Simulate_Time_Step
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Simulate_Time_Step(FLUID_BC<TV>* bc,T time,T dt)
{
    driver->Advance_One_Time_Step(dt);
}

//#####################################################################
// Function Predict_Time_Step
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Predict_Time_Step(T time,T dt)
{
}

//#####################################################################
// Function Before_Time_Step
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Before_Time_Step(T time)
{
    driver->Setup_Fluids(time);
}
//#####################################################################
// Function After_Time_Step
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
After_Time_Step(T time,T dt)
{
}
//#####################################################################
// Function Before_Frame
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Before_Frame(int frame)
{
    driver->Preprocess_Frame(frame+1);
}
//#####################################################################
// Function After_Frame
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
After_Frame(int frame)
{
    driver->Postprocess_Frame(frame+1);
}
//#####################################################################
// Function Make_State
//#####################################################################
template<class TV> FLUID_STATE<TV>* FLUID_SOLVER_PB<TV>::
Make_State() const
{
    return new FLUID_STATE_PB<TV>;
}

//#####################################################################
// Function Make_BC
//#####################################################################
template<class TV> FLUID_BC<TV>* FLUID_SOLVER_PB<TV>::
Make_BC() const
{
    return new FLUID_BC_PB<TV>;
}

//#####################################################################
// Function Save
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Save(FLUID_STATE<TV>* fluid_state) const
{
    FLUID_STATE_PB<TV>& st=dynamic_cast<FLUID_STATE_PB<TV>&>(*fluid_state);
}

//#####################################################################
// Function Restore
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Restore(const FLUID_STATE<TV>* fluid_state)
{
    const FLUID_STATE_PB<TV>& st=dynamic_cast<const FLUID_STATE_PB<TV>&>(*fluid_state);
}

//#####################################################################
// Function Diff_u
//#####################################################################
template<class TV> auto FLUID_SOLVER_PB<TV>::
Diff_u(const FLUID_STATE<TV>* fluid_state) const -> T
{
    return 0;
}

//#####################################################################
// Function Diff_p
//#####################################################################
template<class TV> auto FLUID_SOLVER_PB<TV>::
Diff_p(const FLUID_STATE<TV>* fluid_state) const -> T
{
    return 0;
}

template class FLUID_SOLVER_PB<VECTOR<float,2> >;
template class FLUID_SOLVER_PB<VECTOR<float,3> >;
template class FLUID_SOLVER_PB<VECTOR<double,2> >;
template class FLUID_SOLVER_PB<VECTOR<double,3> >;
}
