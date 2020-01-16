//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Vectors/TWIST.h>
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "SOLID_BC_PB.h"
#include "SOLID_BOUNDARY_VECTOR_PB.h"
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
    driver->Write_Output_Files();
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
Simulate_Time_Step(SOLID_BOUNDARY_VECTOR<TV>* force,T time,T dt)
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

template<class T> T Max_Rot_Disp(const ROTATION<VECTOR<T,3> >& R){return 2*R.Quaternion().v.Magnitude();}
template<class T> T Max_Rot_Disp(const ROTATION<VECTOR<T,2> >& R){return abs(R.Complex()-(T)1);}
template<class T> T Max_Rot_Disp(const ROTATION<VECTOR<T,1> >& R){return 0;}

//#####################################################################
// Function Diff_x
//#####################################################################
template<class TV> auto SOLID_SOLVER_PB<TV>::
Diff_x(const SOLID_STATE<TV>* solid_state) const -> T
{
    const SOLID_STATE_PB<TV>& st=dynamic_cast<const SOLID_STATE_PB<TV>&>(*solid_state);
    T diff_d=(driver->example.solid_body_collection.deformable_body_collection.particles.X-st.X).Maximum_Magnitude();
    auto& frame=driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.frame;
    auto& rb=driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.rigid_body;

    T diff_r=0;
    for(int i=0;i<st.frame.m;i++)
    {
        const FRAME<TV> &f1=frame(i),&f0=st.frame(i);
        T dist_t=(f1.t-f0.t).Magnitude();
        T dist_r=rb(i)->Object_Space_Bounding_Box().Edge_Lengths().Magnitude()*Max_Rot_Disp(f1.r*f0.r.Inverse());
        diff_r=max(diff_r,dist_t+dist_r);
    }

    return max(diff_d,diff_r);
}

//#####################################################################
// Function Diff_v
//#####################################################################
template<class TV> auto SOLID_SOLVER_PB<TV>::
Diff_v(const SOLID_STATE<TV>* solid_state) const -> T
{
    const SOLID_STATE_PB<TV>& st=dynamic_cast<const SOLID_STATE_PB<TV>&>(*solid_state);
    T diff_d=(driver->example.solid_body_collection.deformable_body_collection.particles.V-st.V).Maximum_Magnitude();
    auto& twist=driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.twist;
    auto& rb=driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.rigid_body;

    T diff_r=0;
    for(int i=0;i<st.twist.m;i++)
    {
        const TWIST<TV> &f1=twist(i),&f0=st.twist(i);
        T dist_t=(f1.linear-f0.linear).Magnitude();
        T dist_r=rb(i)->Object_Space_Bounding_Box().Edge_Lengths().Magnitude()*(f1.angular-f0.angular).Magnitude();
        diff_r=max(diff_r,dist_t+dist_r);
    }

    return max(diff_d,diff_r);
}

//#####################################################################
// Function Make_Boundary_Vector
//#####################################################################
template<class TV> SOLID_BOUNDARY_VECTOR<TV>* SOLID_SOLVER_PB<TV>::
Make_Boundary_Vector() const
{
    return new SOLID_BOUNDARY_VECTOR_PB<TV>;
}

//#####################################################################
// Function Get_Velocity
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Get_Velocity(SOLID_BOUNDARY_VECTOR<TV>* v) const
{
    SOLID_BOUNDARY_VECTOR_PB<TV>& bv=dynamic_cast<SOLID_BOUNDARY_VECTOR_PB<TV>&>(*v);
    const auto& V=driver->example.solid_body_collection.deformable_body_collection.particles.V;
    const auto& twist=driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.twist;

    for(auto& p:bv.V)
        p.data=V(p.key);

    for(auto& p:bv.twist)
        p.data=twist(p.key);
}

//#####################################################################
// Function Apply_Velocity_Change
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Apply_Velocity_Change(T c,const SOLID_BOUNDARY_VECTOR<TV>* v) const
{
    const SOLID_BOUNDARY_VECTOR_PB<TV>& bv=dynamic_cast<const SOLID_BOUNDARY_VECTOR_PB<TV>&>(*v);
    auto& V=driver->example.solid_body_collection.deformable_body_collection.particles.V;
    auto& twist=driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.twist;

    for(const auto& p:bv.V)
        V(p.key)+=c*p.data;

    for(const auto& p:bv.twist)
        twist(p.key)+=c*p.data;
}

//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> auto SOLID_SOLVER_PB<TV>::
Inner_Product(const SOLID_BOUNDARY_VECTOR<TV>* u, const SOLID_BOUNDARY_VECTOR<TV>* v) -> T
{
    const SOLID_BOUNDARY_VECTOR_PB<TV>& bu=dynamic_cast<const SOLID_BOUNDARY_VECTOR_PB<TV>&>(*u);
    const SOLID_BOUNDARY_VECTOR_PB<TV>& bv=dynamic_cast<const SOLID_BOUNDARY_VECTOR_PB<TV>&>(*v);
    T x=0;
    TV y;
    TWIST<TV> z;

    for(const auto& p:bv.V)
        if(bu.V.Get(p.key,y))
            x+=p.data.Dot(y);

    for(const auto& p:bv.twist)
        if(bu.twist.Get(p.key,z))
            x+=p.data.Dot(z);

    return x;
}

//#####################################################################
// Function Mass_Inverse
//#####################################################################
template<class TV> void SOLID_SOLVER_PB<TV>::
Mass_Inverse(SOLID_BOUNDARY_VECTOR<TV>* u, const SOLID_BOUNDARY_VECTOR<TV>* v)
{
    SOLID_BOUNDARY_VECTOR_PB<TV>& bu=dynamic_cast<SOLID_BOUNDARY_VECTOR_PB<TV>&>(*u);
    const SOLID_BOUNDARY_VECTOR_PB<TV>& bv=dynamic_cast<const SOLID_BOUNDARY_VECTOR_PB<TV>&>(*v);
    bu.V.Remove_All();
    bu.twist.Remove_All();

    const auto& mi=driver->example.solid_body_collection.deformable_body_collection.particles.one_over_mass;
    const auto& rb=driver->example.solid_body_collection.rigid_body_collection.rigid_body_particles.rigid_body;

    for(const auto& p:bv.V)
        bu.V.Set(p.key,mi(p.key)*p.data);

    for(const auto& p:bv.twist)
        bu.twist.Set(p.key,rb(p.key)->Inertia_Inverse_Times(p.data));
}

template class SOLID_SOLVER_PB<VECTOR<float,2> >;
template class SOLID_SOLVER_PB<VECTOR<float,3> >;
template class SOLID_SOLVER_PB<VECTOR<double,2> >;
template class SOLID_SOLVER_PB<VECTOR<double,3> >;
}
