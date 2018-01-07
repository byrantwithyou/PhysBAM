//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/MATRIX.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Forces_And_Torques/RIGID_PENALTY_WITH_FRICTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_PENALTY_WITH_FRICTION<TV>::
RIGID_PENALTY_WITH_FRICTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
    T stiffness_coefficient,T friction)
    :BASE(rigid_body_collection_input),
    stiffness_coefficient(stiffness_coefficient),friction(friction)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_PENALTY_WITH_FRICTION<TV>::
~RIGID_PENALTY_WITH_FRICTION()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rb1=rigid_body_collection.Rigid_Body(c.b1),
            &rb2=rigid_body_collection.Rigid_Body(c.b2);
        if(c.active){
            TV X=rb1.Frame()*rb1.simplicial_object->particles.X(c.v);
            TV j=stiffness_coefficient*(X-c.Y);
            rigid_F(c.b1)-=rb1.Gather(TWIST<TV>(j,typename TV::SPIN()),X);
            rigid_F(c.b2)+=rb2.Gather(TWIST<TV>(j,typename TV::SPIN()),c.Y);}}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Update_Position_Based_State(const T time)
{
    get_candidates();

    for(int i=0;i<collision_pairs.m;i++)
        Relax_Attachment(i);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,
    ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rb1=rigid_body_collection.Rigid_Body(c.b1),
            &rb2=rigid_body_collection.Rigid_Body(c.b2);
        if(c.active){
            TV X=rb1.Frame()*rb1.simplicial_object->particles.X(c.v);
            TV j=stiffness_coefficient*(rb1.Scatter(rigid_V(c.b1),X).linear-rb2.Scatter(rigid_V(c.b2),c.Y).linear);
            rigid_F(c.b1)-=rb1.Gather(TWIST<TV>(j,typename TV::SPIN()),X);
            rigid_F(c.b2)+=rb2.Gather(TWIST<TV>(j,typename TV::SPIN()),c.Y);}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_PENALTY_WITH_FRICTION<TV>::
Potential_Energy(const T time) const
{
    T pe=0;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active){
            const RIGID_BODY<TV>& rb1=rigid_body_collection.Rigid_Body(c.b1);
            TV X=rb1.Frame()*rb1.simplicial_object->particles.X(c.v);
            pe+=(T).5*stiffness_coefficient*(X-c.Y).Magnitude_Squared();}}
    return pe;
}
//#####################################################################
// Function Relax_Attachment
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Relax_Attachment(int cp)
{
    COLLISION_PAIR& c=collision_pairs(cp);
    const RIGID_BODY<TV>& rb1=rigid_body_collection.Rigid_Body(c.b1),
        &rb2=rigid_body_collection.Rigid_Body(c.b2);
    const IMPLICIT_OBJECT<TV>* io=rb2.implicit_object;
    TV X=rb2.Frame()*c.X;
    TV Z=rb1.Frame()*rb1.simplicial_object->particles.X(c.v);
    T phi=io->Extended_Phi(Z);
    TV n=io->Extended_Normal(Z);
    auto pr=Relax_Attachment_Helper(Z,X,phi,n,friction);
    c.Y=pr.Y;
    c.active=pr.active;
    // TODO: Fix derivatives
}
//#####################################################################
// Function Update_Attachments_And_Prune_Pairs
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Update_Attachments_And_Prune_Pairs()
{
    int k=0;
    for(int i=0;i<collision_pairs.m;i++){
        COLLISION_PAIR c=collision_pairs(i);
        if(c.active){
            const RIGID_BODY<TV>& rb2=rigid_body_collection.Rigid_Body(c.b2);
            c.X=rb2.Frame().Inverse_Times(c.Y);
            collision_pairs(k++)=c;}
        else hash.Delete({{c.b1,c.v},c.b2});}
    collision_pairs.Resize(k);
}
//#####################################################################
// Function Add_Pair
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Add_Pair(int b1,int v,int b2)
{
    // TODO: Interpolate X^n and X^(n+1) to choose surface point.
    if(hash.Contains({{b1,v},b2})) return;
    const RIGID_BODY<TV>& rb1=rigid_body_collection.Rigid_Body(b1),
        &rb2=rigid_body_collection.Rigid_Body(b2);
    TV X=rb1.Frame()*rb1.simplicial_object->particles.X(v);
    if(rb2.implicit_object->Extended_Phi(X)>0) return;
    TV W=rb2.implicit_object->Closest_Point_On_Boundary(X);
    COLLISION_PAIR c={b1,v,b2,rb2.Frame().Inverse_Times(W)};
    collision_pairs.Append(c);
    hash.Insert({{b1,v},b2});
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,
    ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated)
{
}
//#####################################################################
// Function Use_Rest_State_For_Strain_Rate
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Limit_Time_Step_By_Strain_Rate
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input,const T max_strain_per_time_step_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> rigid_frequency)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int RIGID_PENALTY_WITH_FRICTION<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_PENALTY_WITH_FRICTION<TV>::
CFL_Strain_Rate() const
{
    PHYSBAM_FATAL_ERROR();
}
namespace PhysBAM{
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<float,2> >;
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<double,2> >;
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<float,3> >;
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<double,3> >;
}
