//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/COMBINED_COLLISIONS_COUPLED_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
COMBINED_COLLISIONS_COUPLED_IMPULSE(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions,bool do_collision)
    :DBASE(solid_body_collection_input.deformable_body_collection),RBASE(solid_body_collection_input.rigid_body_collection,rigid_body_collisions,do_collision)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
~COMBINED_COLLISIONS_COUPLED_IMPULSE()
{
}
//#####################################################################
// Function Resize
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
Resize()
{
    DBASE::Resize();
    RBASE::Resize();
}
//#####################################################################
// Function Apply
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
Apply(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)
{
    for(int i=1;i<=list.m;i++){
        if(Combined_Body_Id_Is_Rigid(list(i)))
            RBASE::Apply(Combined_Body_Id_To_Rigid_Body(list(i)),dt,time);
        else
            DBASE::Apply(Combined_Body_Id_To_Deformable_Particle(list(i)),dt,time);}
}
//#####################################################################
// Function Clear
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
Clear(const ARRAY<COMBINED_BODY_ID>& list)
{
    for(int i=1;i<=list.m;i++){
        if(Combined_Body_Id_Is_Rigid(list(i))){
            int p=Combined_Body_Id_To_Rigid_Body(list(i));
            wrench(p)=TWIST<TV>();}
        else{
            int p=Combined_Body_Id_To_Deformable_Particle(list(i));
            impulse(p)=TV();}}
}
//#####################################################################
// Function Clear
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
Clear()
{
    DBASE::Clear();
    RBASE::Clear();
}
//#####################################################################
// Function Apply_Coupled_Impulse
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
Apply_Coupled_Impulse(int p,int b,const TV& j)
{
    impulse(p)+=j;
    Apply_Rigid_Impulse(b,DBASE::deformable_body_collection.particles.X(p),TWIST<TV>(-j,typename TV::SPIN()));
}
//#####################################################################
// Function Setup_Discover_State
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
Setup_Discover_State(T dt,T time)
{
    DBASE::Setup_Discover_State(dt,time);
    RBASE::Setup_Discover_State(dt,time);
}
//#####################################################################
// Function Setup_Impulse_State
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
Setup_Impulse_State(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)
{
    for(int i=1;i<=list.m;i++){
        if(Combined_Body_Id_Is_Rigid(list(i))){
            int p=Combined_Body_Id_To_Rigid_Body(list(i));
            RBASE::Setup_Impulse_State(p,dt,time);}
        else{
            int p=Combined_Body_Id_To_Deformable_Particle(list(i));
            DBASE::Setup_Impulse_State(p,dt,time);}}
}
//#####################################################################
// Function Finish_Step
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
Finish_Step(const ARRAY<COMBINED_BODY_ID>& list,T dt,T time)
{
    for(int i=1;i<=list.m;i++){
        if(Combined_Body_Id_Is_Rigid(list(i))){
            int p=Combined_Body_Id_To_Rigid_Body(list(i));
            RBASE::Finish_Step(p,dt,time);}
        else{
            int p=Combined_Body_Id_To_Deformable_Particle(list(i));
            DBASE::Finish_Step(p,dt,time);}}
}
//#####################################################################
// Function Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>::
Kinetic_Energy(const ARRAY<COMBINED_BODY_ID>& list) const
{
    T ke=0;
    for(int i=1;i<=list.m;i++){
        if(Combined_Body_Id_Is_Rigid(list(i)))
            ke+=RBASE::Kinetic_Energy(Combined_Body_Id_To_Rigid_Body(list(i)));
        else
            ke+=DBASE::Kinetic_Energy(Combined_Body_Id_To_Deformable_Particle(list(i)));}
    return ke;
}
template class COMBINED_COLLISIONS_COUPLED_IMPULSE<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_COUPLED_IMPULSE<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_COUPLED_IMPULSE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_COUPLED_IMPULSE<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_COUPLED_IMPULSE<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_COUPLED_IMPULSE<VECTOR<double,3> >;
#endif
