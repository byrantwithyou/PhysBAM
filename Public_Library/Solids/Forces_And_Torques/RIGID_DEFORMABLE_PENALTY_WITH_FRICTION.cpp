//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/FINE_TIMER.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Matrices/MATRIX.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Rigids/Forces_And_Torques/MOVE_RIGID_BODY_DIFF.h>
#include <Rigids/Forces_And_Torques/RIGID_PENALTY_WITH_FRICTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/RIGID_DEFORMABLE_PENALTY_WITH_FRICTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
RIGID_DEFORMABLE_PENALTY_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input,
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
    const ARRAY<MOVE_RIGID_BODY_DIFF<TV> >& move_rb_diff)
    :BASE(particles_input,rigid_body_collection_input),
    move_rb_diff(move_rb_diff)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
~RIGID_DEFORMABLE_PENALTY_WITH_FRICTION()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    TIMER_SCOPE_FUNC;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(c.b);
        if(c.active){
            TV j=stiffness_coefficient*(particles.X(c.p)-c.Y);
            F(c.p)-=j;
            if(!rb.Has_Infinite_Inertia())
                rigid_F(c.b)+=rb.Gather(TWIST<TV>(j,typename TV::SPIN()),c.Y);}}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Update_Position_Based_State(const T time)
{
    TIMER_SCOPE_FUNC;
    get_candidates();

    num_dynamic=0;
    num_stick=0;
    for(int i=0;i<collision_pairs.m;i++)
        Relax_Attachment(i);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,
    ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,
    ARRAY_VIEW<TWIST<TV> > rigid_F,const T time,bool transpose) const
{
    TIMER_SCOPE_FUNC;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(c.b);
        if(c.active){
            TV j=stiffness_coefficient*(particles.X(c.p)-c.Y);
            if(rb.Has_Infinite_Inertia()){
                if(transpose)
                    F(c.p)+=(c.dYdZ-1).Transpose_Times(stiffness_coefficient*V(c.p));
                else F(c.p)+=(c.dYdZ-1)*(stiffness_coefficient*V(c.p));}
            else{
                if(transpose){
                    TWIST<TV> tw=rb.Scatter(rigid_V(c.b),c.Y);
                    TV dj=tw.linear-V(c.p),cp=j.Cross(tw.angular);
                    TV dZ=stiffness_coefficient*dj,dY=cp-dZ;
                    rigid_F(c.b).linear+=c.dYdL.Transpose_Times(dY)-cp;
                    rigid_F(c.b).angular+=c.dYdA.Transpose_Times(dY);
                    F(c.p)+=dZ+c.dYdZ.Transpose_Times(dY);}
                else{
                    TV dZ=V(c.p),dL=rigid_V(c.b).linear;
                    auto dA=rigid_V(c.b).angular;
                    TV dY=c.dYdZ*dZ+c.dYdL*dL+c.dYdA*dA;
                    TV dj=stiffness_coefficient*(dZ-dY);
                    F(c.p)-=dj;
                    rigid_F(c.b)+=rb.Gather(TWIST<TV>(dj,(dY-dL).Cross(j)),c.Y);}}}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Potential_Energy(const T time) const
{
    TIMER_SCOPE_FUNC;
    T pe=0;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active)
            pe+=(T).5*stiffness_coefficient*(particles.X(c.p)-c.Y).Magnitude_Squared();}
    return pe;
}
//#####################################################################
// Function Relax_Attachment
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Relax_Attachment(int cp)
{
    TIMER_SCOPE_FUNC;
    COLLISION_PAIR& c=collision_pairs(cp);
    const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(c.b);
    const IMPLICIT_OBJECT<TV>* io=rb.implicit_object->object_space_implicit_object;
    const MOVE_RIGID_BODY_DIFF<TV>& mr=move_rb_diff(c.b);
    MATRIX<T,TV::m> dXdv,dXdL,dZdv,dWdZ,dWdL,dVdY,dVdL;
    MATRIX<T,TV::m,TV::SPIN::m> dXdA,dWdA,dVdA;

    TV Z=particles.X(c.p);
    TV X=mr.Frame_Times(c.X,dXdv,dXdL,dXdA),W,V;

    c.active=Project_Attachment_To_Surface(W,mr,io,Z,dWdZ,dWdL,dWdA,true);
    if(!c.active) return;

    auto pr=Relax_Attachment_Helper(Z,X,W,friction);
    if(pr.dynamic) num_dynamic++;
    else num_stick++;
    
    Project_Attachment_To_Surface(V,mr,io,pr.Y,dVdY,dVdL,dVdA,false);
    c.Y=V;
    MATRIX<T,TV::m> dVdZ=dVdY*(pr.dYdZ+pr.dYdW*dWdZ);
    c.dYdZ=dVdZ;
    c.dYdL=dVdY*(pr.dYdX*dXdL+pr.dYdW*dWdL)+dVdL;
    c.dYdA=dVdY*(pr.dYdX*dXdA+pr.dYdW*dWdA)+dVdA;
}
//#####################################################################
// Function Update_Attachments_And_Prune_Pairs
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Update_Attachments_And_Prune_Pairs()
{
    TIMER_SCOPE_FUNC;
    int k=0;
    for(int i=0;i<collision_pairs.m;i++){
        COLLISION_PAIR c=collision_pairs(i);
        if(c.active){
            const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(c.b);
            c.X=rb.Frame().Inverse_Times(c.Y);
            collision_pairs(k++)=c;}
        else hash.Delete({c.p,c.b});}
    collision_pairs.Resize(k);
}
//#####################################################################
// Function Add_Pair
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Pair(int p,int b)
{
    TIMER_SCOPE_FUNC;
    // TODO: Interpolate X^n and X^(n+1) to choose surface point.
    if(hash.Contains({p,b})) return;
    TV X=particles.X(p);
    const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(b);
    if(rb.implicit_object->Extended_Phi(X)>0) return;
    TV W=rb.implicit_object->Closest_Point_On_Boundary(X);
    COLLISION_PAIR c={p,b,rb.Frame().Inverse_Times(W)};
    collision_pairs.Append(c);
    hash.Insert({p,b});
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,
    ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,
    ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,
    const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Use_Rest_State_For_Strain_Rate
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Limit_Time_Step_By_Strain_Rate
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input,const T max_strain_per_time_step_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::DEFORMABLE_FREQUENCY_DATA> frequency,ARRAY_VIEW<typename BASE::RIGID_FREQUENCY_DATA> rigid_frequency)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Force_Data
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
CFL_Strain_Rate() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Read(TYPED_ISTREAM& input)
{
    Read_Binary(input,collision_pairs);
    ARRAY<PAIR<int,int> > keys;
    Read_Binary(input,keys);
    hash.Set_All(keys);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,collision_pairs);
    ARRAY<PAIR<int,int> > keys;
    hash.Get_Keys(keys);
    keys.Sort();
    Write_Binary(output,keys);
}
namespace PhysBAM{
template class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<VECTOR<float,2> >;
template class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<VECTOR<double,2> >;
template class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<VECTOR<float,3> >;
template class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<VECTOR<double,3> >;
}
