//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/COMBINED_COLLISIONS_RIGID_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
COMBINED_COLLISIONS_RIGID_IMPULSE(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions_input,bool do_collision)
    :rigid_body_collection(rigid_body_collection_input),rigid_body_collisions(rigid_body_collisions_input),update_positions_for_collisions_on_apply(do_collision),
    update_positions_for_contact_on_apply(!do_collision),use_collisions(do_collision)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
~COMBINED_COLLISIONS_RIGID_IMPULSE()
{
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Resize()
{
    wrench.Resize(rigid_body_collection.rigid_geometry_collection.particles.array_collection->Size());
}
//#####################################################################
// Function Apply
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Apply(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)
{
    for(int i=1;i<=list.m;i++)
        Apply(Combined_Body_Id_To_Rigid_Body(list(i)),dt,time);
}
//#####################################################################
// Function Apply
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Apply(int b,T dt,T time)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(b);
    if(rigid_body.Has_Infinite_Inertia()) return;
    rigid_body.Apply_Impulse_To_Body(rigid_body.X(),wrench(b).linear,wrench(b).angular);
}
//#####################################################################
// Function Clear
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Clear(const ARRAY<COMBINED_BODY_ID>& list)
{
    for(int i=1;i<=list.m;i++){
        int p=Combined_Body_Id_To_Rigid_Body(list(i));
        wrench(p)=TWIST<TV>();}
}
//#####################################################################
// Function Clear
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Clear()
{
    ARRAYS_COMPUTATIONS::Fill(wrench,TWIST<TV>());
}
//#####################################################################
// Function Apply_Rigid_Impulse
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Apply_Rigid_Impulse(int b,const TV& location,const TWIST<TV>& j)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(b);
    if(rigid_body.Has_Infinite_Inertia()) return;
    wrench(b)+=rigid_body.Gather(j,location);
}
//#####################################################################
// Function Apply_Rigid_Impulse
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Apply_Rigid_Impulse(int a,int b,const TV& location,const TWIST<TV>& j)
{
    Apply_Rigid_Impulse(a,location,j);
    Apply_Rigid_Impulse(b,location,-j);
}
//#####################################################################
// Function Body_Twist
//#####################################################################
template<class TV> TWIST<TV> COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Body_Twist(int b) const
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(b);
    return rigid_body.Twist()+rigid_body.Inertia_Inverse_Times(wrench(b));
}
//#####################################################################
// Function Velocity_At_Point
//#####################################################################
template<class TV> TV COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Velocity_At_Point(int e,const TV& location) const
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(e);
    if(rigid_body.Has_Infinite_Inertia()) return rigid_body.Pointwise_Object_Velocity(location);
    return rigid_body.Pointwise_Object_Velocity(Body_Twist(e),rigid_body.X(),location);
}
//#####################################################################
// Function Setup_Discover_State
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Setup_Discover_State(T dt,T time)
{
}
//#####################################################################
// Function Setup_Impulse_State
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Setup_Impulse_State(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)
{
    if(use_collisions) return;
    for(int i=1;i<=list.m;i++)
        Setup_Impulse_State(Combined_Body_Id_To_Rigid_Body(list(i)),dt,time);
}
//#####################################################################
// Function Finish_Step
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Finish_Step(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)
{
    if(use_collisions) return;
    for(int i=1;i<=list.m;i++)
        Finish_Step(Combined_Body_Id_To_Rigid_Body(list(i)),dt,time);
}
//#####################################################################
// Function Setup_Impulse_State
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Setup_Impulse_State(int b,T dt,T time)
{
    if(!use_collisions)
        rigid_body_collisions.collision_callbacks.Swap_State(b);
}
//#####################################################################
// Function Finish_Step
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Finish_Step(int b,T dt,T time)
{
    if(!use_collisions){
        rigid_body_collisions.collision_callbacks.Save_Position(b);
        rigid_body_collisions.Euler_Step_Position(b,dt,time);}
}
//#####################################################################
// Function Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Kinetic_Energy(const ARRAY<COMBINED_BODY_ID>& list) const
{
    T ke=0;
    for(int i=1;i<=list.m;i++)
        ke+=Kinetic_Energy(Combined_Body_Id_To_Rigid_Body(list(i)));
    return ke;
}
//#####################################################################
// Function Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Kinetic_Energy(int b) const
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(b);
    if(rigid_body.Has_Infinite_Inertia()) return 0;
    return rigid_body.Kinetic_Energy(Body_Twist(b));
}
//#####################################################################
// Function Kinetic_Energy_Change
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_RIGID_IMPULSE<TV>::
Kinetic_Energy_Change(int b,const TV& location,const TWIST<TV>& j) const
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(b);
    if(rigid_body.Has_Infinite_Inertia()) return 0;
    TWIST<TV> twist=Body_Twist(b),rj=rigid_body.Gather(j,location);
    return TV::Dot_Product(twist.linear,rj.linear)+TV::SPIN::Dot_Product(twist.angular,rj.angular);
}
template class COMBINED_COLLISIONS_RIGID_IMPULSE<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_RIGID_IMPULSE<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_RIGID_IMPULSE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_RIGID_IMPULSE<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_RIGID_IMPULSE<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_RIGID_IMPULSE<VECTOR<double,3> >;
#endif
