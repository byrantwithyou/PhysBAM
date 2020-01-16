//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Projection/BOUNDARY_CONDITION_DOUBLE_FINE.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Fluids/Fluids/FLUID_COLLECTION.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "FLUID_BC_PB.h"
#include "FLUID_BOUNDARY_VECTOR_PB.h"
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
    auto old_cb=driver->example.get_unified_boundary_conditions;
}

//#####################################################################
// Function Write
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Write(int frame) const
{
    driver->Write_Output_Files();
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
Simulate_Time_Step(FLUID_BOUNDARY_VECTOR<TV>* velocity,T time,T dt)
{
    auto* bc_v=dynamic_cast<FLUID_BOUNDARY_VECTOR_PB<TV>*>(velocity);
    auto* bc_fine=driver->example.fluids_parameters.bc_fine;

    ARRAY<char,TV_INT> bc_type,bc_type_current;
    HASHTABLE<TV_INT,PAIR<T,int> > bc_p;
    HASHTABLE<FACE_INDEX<TV::m>,PAIR<T,int> > bc_u;

    T bc_v_value=0;
    for(auto& h:bc_fine->bc_u)
        if(bc_v->V.Get(h.key,bc_v_value))
            h.data={bc_v_value,1};
            
    bc_type.Exchange(bc_fine->bc_type);
    bc_type_current.Exchange(bc_fine->bc_type_current);
    bc_p.Exchange(bc_fine->bc_p);
    bc_u.Exchange(bc_fine->bc_u);
    auto old_cb=driver->example.get_unified_boundary_conditions;
    driver->example.get_unified_boundary_conditions=
        [&,this](BOUNDARY_CONDITION_DOUBLE_FINE<TV>* bc_fine,const T time)
        {
            bc_type.Exchange(bc_fine->bc_type);
            bc_type_current.Exchange(bc_fine->bc_type_current);
            bc_p.Exchange(bc_fine->bc_p);
            bc_u.Exchange(bc_fine->bc_u);
        };
    driver->Advance_One_Time_Step(dt);
    driver->example.get_unified_boundary_conditions=old_cb;
}

//#####################################################################
// Function Predict_Time_Step
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Predict_Time_Step(T time,T dt)
{
    // Predict that the state does not change.
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
    st.face_velocities=driver->example.fluid_collection.incompressible_fluid_collection.face_velocities;
    st.pressure=driver->example.fluids_parameters.incompressible->projection.p;
}

//#####################################################################
// Function Restore
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Restore(const FLUID_STATE<TV>* fluid_state)
{
    const FLUID_STATE_PB<TV>& st=dynamic_cast<const FLUID_STATE_PB<TV>&>(*fluid_state);
    driver->example.fluid_collection.incompressible_fluid_collection.face_velocities=st.face_velocities;
    driver->example.fluids_parameters.incompressible->projection.p=st.pressure;
}

//#####################################################################
// Function Diff_u
//#####################################################################
template<class TV> auto FLUID_SOLVER_PB<TV>::
Diff_u(const FLUID_STATE<TV>* fluid_state) const -> T
{
    const FLUID_STATE_PB<TV>& st=dynamic_cast<const FLUID_STATE_PB<TV>&>(*fluid_state);
    return (driver->example.fluid_collection.incompressible_fluid_collection.face_velocities.array-st.face_velocities.array).Max_Abs();
}

//#####################################################################
// Function Diff_p
//#####################################################################
template<class TV> auto FLUID_SOLVER_PB<TV>::
Diff_p(const FLUID_STATE<TV>* fluid_state) const -> T
{
    const FLUID_STATE_PB<TV>& st=dynamic_cast<const FLUID_STATE_PB<TV>&>(*fluid_state);
    return (driver->example.fluids_parameters.incompressible->projection.p.array-st.pressure.array).Max_Abs();
}

//#####################################################################
// Function Get_Constraints
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Get_Constraints(ARRAY<FLUID_BOUNDARY_VECTOR<TV>*>& array) const
{
}

//#####################################################################
// Function Make_Boundary_Vector
//#####################################################################
template<class TV> FLUID_BOUNDARY_VECTOR<TV>* FLUID_SOLVER_PB<TV>::
Make_Boundary_Vector() const
{
    return new FLUID_BOUNDARY_VECTOR_PB<TV>;
}

//#####################################################################
// Function Get_Force
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Get_Force(FLUID_BOUNDARY_VECTOR<TV>* force) const
{
}

template class FLUID_SOLVER_PB<VECTOR<float,2> >;
template class FLUID_SOLVER_PB<VECTOR<float,3> >;
template class FLUID_SOLVER_PB<VECTOR<double,2> >;
template class FLUID_SOLVER_PB<VECTOR<double,3> >;
}
