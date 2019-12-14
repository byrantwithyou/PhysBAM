//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "SOLID_BC_PB.h"
#include "SOLID_SOLVER_PB.h"
#include "SOLID_STATE_PB.h"
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLID_SOLVER_PB<TV>::
SOLID_SOLVER_PB()
{
}

//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLID_SOLVER_PB<TV>::
~SOLID_SOLVER_PB()
{
}

//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Initialize()
{
    driver->Initialize();
}

//#####################################################################
// Function Write
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Write(int frame) const
{
    driver->Write_Output_Files(frame);
}

//#####################################################################
// Function Read
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Read(int frame)
{
}

//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> auto SOLID_SOLVER_PB<TV>::
Compute_Dt(T time) const -> T
{
    return driver->Compute_Solids_Dt(time);
}

//#####################################################################
// Function Simulate_Time_Step
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Simulate_Time_Step(SOLID_BC<TV>* bc,T time,T dt)
{
    driver->Advance_One_Time_Step(dt);
}
//#####################################################################
// Function Before_Time_Step
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Before_Time_Step(T time)
{
    driver->Setup_Solids(time);
}
//#####################################################################
// Function After_Time_Step
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
After_Time_Step(T time,T dt)
{
}
//#####################################################################
// Function Before_Frame
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Before_Frame(int frame)
{
    driver->Preprocess_Frame(frame+1);
}
//#####################################################################
// Function After_Frame
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
After_Frame(int frame)
{
    driver->Postprocess_Frame(frame+1);
}
//#####################################################################
// Function Make_State
//#####################################################################
template<class TV> SOLID_STATE<TV>* SOLID_SOLVER_PB<TV>::
Make_State() const
{
    return new SOLID_STATE_PB<TV>;
}

//#####################################################################
// Function Make_BC
//#####################################################################
template<class TV> SOLID_BC<TV>* SOLID_SOLVER_PB<TV>::
Make_BC() const
{
    return new SOLID_BC_PB<TV>;
}

//#####################################################################
// Function Save
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Save(SOLID_STATE<TV>* solid_state) const
{
    SOLID_STATE_PB<TV>& st=dynamic_cast<SOLID_STATE_PB<TV>&>(*solid_state);
    st.time=driver->time;
    st.frame=driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.frame;
    st.twist=driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.twist;
    st.X=driver->example.solid_body_collection.deformable_body_collection.particles.X;
    st.V=driver->example.solid_body_collection.deformable_body_collection.particles.V;
    st.repulsion_pair_update_count=
        driver->example.solids_parameters.triangle_collision_parameters.repulsion_pair_update_count;
    st.topological_hierarchy_build_count=
        driver->example.solids_parameters.triangle_collision_parameters.topological_hierarchy_build_count;
}

//#####################################################################
// Function Restore
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Restore(const SOLID_STATE<TV>* solid_state)
{
    const SOLID_STATE_PB<TV>& st=dynamic_cast<const SOLID_STATE_PB<TV>&>(*solid_state);
    driver->time=st.time;
    driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.frame=st.frame;
    driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.twist=st.twist;
    driver->example.solid_body_collection.rigid_body_collection.Update_Angular_Momentum();
    driver->example.solid_body_collection.deformable_body_collection.particles.X=st.X;
    driver->example.solid_body_collection.deformable_body_collection.particles.V=st.V;
    driver->example.solids_parameters.triangle_collision_parameters.repulsion_pair_update_count=
        st.repulsion_pair_update_count;
    driver->example.solids_parameters.triangle_collision_parameters.topological_hierarchy_build_count=
        st.topological_hierarchy_build_count;
}

//#####################################################################
// Function Diff_x
//#####################################################################
template<class TV> auto SOLID_SOLVER_PB<TV>::
Diff_x(const SOLID_STATE<TV>* solid_state) const -> T
{
    return 0;
}

//#####################################################################
// Function Diff_v
//#####################################################################
template<class TV> auto SOLID_SOLVER_PB<TV>::
Diff_v(const SOLID_STATE<TV>* solid_state) const -> T
{
    return 0;
}

template class SOLID_SOLVER_PB<VECTOR<float,2> >;
template class SOLID_SOLVER_PB<VECTOR<float,3> >;
template class SOLID_SOLVER_PB<VECTOR<double,2> >;
template class SOLID_SOLVER_PB<VECTOR<double,3> >;
}
