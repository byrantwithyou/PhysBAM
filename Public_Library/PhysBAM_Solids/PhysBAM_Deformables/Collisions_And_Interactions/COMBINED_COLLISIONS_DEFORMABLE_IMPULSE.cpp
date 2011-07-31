//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_DEFORMABLE_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
COMBINED_COLLISIONS_DEFORMABLE_IMPULSE(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input)
    :deformable_body_collection(deformable_body_collection_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
~COMBINED_COLLISIONS_DEFORMABLE_IMPULSE()
{
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Resize()
{
    impulse.Resize(deformable_body_collection.particles.array_collection->Size());
}
//#####################################################################
// Function Apply
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Apply(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)
{
    for(int i=1;i<=list.m;i++)
        Apply(Combined_Body_Id_To_Deformable_Particle(list(i)),dt,time);
}
//#####################################################################
// Function Apply
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Apply(int e,T dt,T time)
{
    TV dv=impulse(e)*deformable_body_collection.particles.one_over_mass(e);
    deformable_body_collection.particles.V(e)+=dv;
    deformable_body_collection.particles.X(e)+=dv*dt;
}
//#####################################################################
// Function Clear
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Clear(const ARRAY<COMBINED_BODY_ID>& list)
{
    for(int i=1;i<=list.m;i++){
        int p=Combined_Body_Id_To_Deformable_Particle(list(i));
        impulse(p)=TV();}
}
//#####################################################################
// Function Clear
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Clear()
{
    ARRAYS_COMPUTATIONS::Fill(impulse,TV());
}
//#####################################################################
// Function Particle_Velocity
//#####################################################################
template<class TV> TV COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Particle_Velocity(int e) const
{
    return deformable_body_collection.particles.V(e)+impulse(e)*deformable_body_collection.particles.one_over_mass(e);
}
//#####################################################################
// Function Setup_Discover_State
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Setup_Discover_State(T dt,T time)
{
}
//#####################################################################
// Function Setup_Impulse_State
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Setup_Impulse_State(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)
{
}
//#####################################################################
// Function Finish_Step
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Finish_Step(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)
{
}
//#####################################################################
// Function Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Kinetic_Energy(const ARRAY<COMBINED_BODY_ID>& list) const
{
    T ke=0;
    for(int i=1;i<=list.m;i++)
        ke+=Kinetic_Energy(Combined_Body_Id_To_Deformable_Particle(list(i)));
    return ke;
}
//#####################################################################
// Function Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Kinetic_Energy(int p) const
{
    return (T).5*deformable_body_collection.particles.mass(p)*Particle_Velocity(p).Magnitude_Squared();
}
//#####################################################################
// Function Kinetic_Energy_Change
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>::
Kinetic_Energy_Change(int b,const TV& j) const
{
    return TV::Dot_Product(Particle_Velocity(b),j);
}
template struct COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<VECTOR<float,1> >;
template struct COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<VECTOR<float,2> >;
template struct COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template struct COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<VECTOR<double,1> >;
template struct COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<VECTOR<double,2> >;
template struct COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<VECTOR<double,3> >;
#endif
