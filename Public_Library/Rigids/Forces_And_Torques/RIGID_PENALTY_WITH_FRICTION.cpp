//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/FINE_TIMER.h>
#include <Core/Matrices/MATRIX.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Forces_And_Torques/MOVE_RIGID_BODY_DIFF.h>
#include <Rigids/Forces_And_Torques/RIGID_PENALTY_WITH_FRICTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_PENALTY_WITH_FRICTION<TV>::
RIGID_PENALTY_WITH_FRICTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
    const ARRAY<MOVE_RIGID_BODY_DIFF<TV> >& move_rb_diff)
    :BASE(rigid_body_collection_input),
    move_rb_diff(move_rb_diff)
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
    TIMER_SCOPE_FUNC;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(c.bs),
            &rbi=rigid_body_collection.Rigid_Body(c.bi);
        bool ts=!rbs.Has_Infinite_Inertia(),ti=!rbi.Has_Infinite_Inertia();
        if(c.active){
            TV j=stiffness_coefficient*(c.Z-c.Y);
            if(ts) rigid_F(c.bs)-=rbs.Gather(TWIST<TV>(j,typename TV::SPIN()),c.Z);
            if(ti) rigid_F(c.bi)+=rbi.Gather(TWIST<TV>(j,typename TV::SPIN()),c.Y);}}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Update_Position_Based_State(const T time)
{
    TIMER_SCOPE_FUNC;
    get_candidates();

    for(int i=0;i<collision_pairs.m;i++)
        Relax_Attachment(i);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,
    ARRAY_VIEW<TWIST<TV> > rigid_F,const T time,bool transpose) const
{
    TIMER_SCOPE_FUNC;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(c.bs),
            &rbi=rigid_body_collection.Rigid_Body(c.bi);
        bool ts=!rbs.Has_Infinite_Inertia(),ti=!rbi.Has_Infinite_Inertia();
        if(c.active){
            TV j=stiffness_coefficient*(c.Z-c.Y);
            if(transpose){
                TWIST<TV> qs,qi;
                if(ts) qs=rbs.Scatter(rigid_V(c.bs),c.Z);
                if(ti) qi=rbi.Scatter(rigid_V(c.bi),c.Y);
                TV dj=stiffness_coefficient*(qi.linear-qs.linear);
                TV cs=j.Cross(qs.angular),ci=j.Cross(qi.angular);
                TV dY=ci-dj,dZ=dj-cs;
                if(ti){
                    rigid_F(c.bi).linear+=c.dYdLi.Transpose_Times(dY)-ci;
                    rigid_F(c.bi).angular+=c.dYdAi.Transpose_Times(dY);}
                if(ts){
                    rigid_F(c.bs).linear+=cs+c.dYdLs.Transpose_Times(dY)+c.dZdLs.Transpose_Times(dZ);
                    rigid_F(c.bs).angular+=c.dYdAs.Transpose_Times(dY)+c.dZdAs.Transpose_Times(dZ);}}
            else{
                TV dLs,dLi;
                typename TV::SPIN dAs,dAi;
                if(ts){
                    dLs=rigid_V(c.bs).linear;
                    dAs=rigid_V(c.bs).angular;}
                if(ti){
                    dLi=rigid_V(c.bi).linear;
                    dAi=rigid_V(c.bi).angular;}
                TV dZ=c.dZdLs*dLs+c.dZdAs*dAs;
                TV dY=c.dYdLs*dLs+c.dYdAs*dAs+c.dYdLi*dLi+c.dYdAi*dAi;
                TV dj=stiffness_coefficient*(dZ-dY);
                if(ts) rigid_F(c.bs)-=rbs.Gather(TWIST<TV>(dj,(dZ-dLs).Cross(j)),c.Z);
                if(ti) rigid_F(c.bi)+=rbi.Gather(TWIST<TV>(dj,(dY-dLi).Cross(j)),c.Y);}}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_PENALTY_WITH_FRICTION<TV>::
Potential_Energy(const T time) const
{
    TIMER_SCOPE_FUNC;
    T pe=0;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active){
            const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(c.bs);
            TV X=rbs.Frame()*rbs.simplicial_object->particles.X(c.v);
            pe+=(T).5*stiffness_coefficient*(X-c.Y).Magnitude_Squared();}}
    return pe;
}
template<class TV> bool MOVING_LEVEL_SET_HELPER<TV>::
Init(const MOVE_RIGID_BODY_DIFF<TV>& mr,const IMPLICIT_OBJECT<TV>* io,const TV& X,bool exit_if_sep)
{
    U=mr.Frame_Inverse_Times(X,dUdX,dUdL,dUdA);
    phi=io->Extended_Phi(U);
    if(exit_if_sep && phi>0) return false;
    N=io->Extended_Normal(U);
    MATRIX<T,TV::m> dndN;
    n=mr.Rotate(N,dndN,dndA);
    SYMMETRIC_MATRIX<T,TV::m> dNdU=io->Hessian(U);
    dndU=dndN*dNdU;
    return true;
}
//#####################################################################
// Function Project_Attachment_To_Surface
//#####################################################################
template<class TV,class T> void PhysBAM::
Project_Attachment_To_Surface(TV& W,MOVING_LEVEL_SET_HELPER<TV>& s,
    const TV& X,MATRIX<T,TV::m>& dWdX,
    MATRIX<T,TV::m>& dWdL,MATRIX<T,TV::m,TV::SPIN::m>& dWdA)
{
    TIMER_SCOPE_FUNC;
    W=X-s.phi*s.n;
    MATRIX<T,TV::m> dWdU=-s.phi*s.dndU-Outer_Product(s.n,s.N);
    dWdX=(T)1+dWdU*s.dUdX;
    dWdL=dWdU*s.dUdL;
    dWdA=dWdU*s.dUdA-s.phi*s.dndA;
}
//#####################################################################
// Function Relax_Attachment
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Relax_Attachment(int cp)
{
    TIMER_SCOPE_FUNC;
    COLLISION_PAIR& c=collision_pairs(cp);
    const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(c.bs),
        &rbi=rigid_body_collection.Rigid_Body(c.bi);
    const IMPLICIT_OBJECT<TV>* io=rbi.implicit_object->object_space_implicit_object;
    const MOVE_RIGID_BODY_DIFF<TV>& mrs=move_rb_diff(c.bs);
    const MOVE_RIGID_BODY_DIFF<TV>& mri=move_rb_diff(c.bi);
    MATRIX<T,TV::m> dXdv,dXdLi,dZdv,dWdZ,dWdLi,dYdK,dYdLi,dKdLi;
    MATRIX<T,TV::m,TV::SPIN::m> dXdAi,dZdAs,dWdAi,dYdAi,dKdAi;

    TV Xs=rbs.simplicial_object->particles.X(c.v);
    c.Z=mrs.Frame_Times(Xs,dZdv,c.dZdLs,c.dZdAs);
    TV X=mri.Frame_Times(c.X,dXdv,dXdLi,dXdAi),W,Y;

    MOVING_LEVEL_SET_HELPER<TV> mZ,mK;
    c.active=mZ.Init(mri,io,c.Z,true);
    if(!c.active) return;
    Project_Attachment_To_Surface(W,mZ,c.Z,dWdZ,dWdLi,dWdAi);

    RELAX_ATTACHMENT_HELPER<TV> h;
    if(use_bisection) Relax_Attachment_Helper_Search(h,c.Z,X,W,rbi.implicit_object,friction);
    else Relax_Attachment_Helper(h,c.Z,X,W,friction);

    mK.Init(mri,io,h.K,false);
    Project_Attachment_To_Surface(Y,mK,h.K,dYdK,dYdLi,dYdAi);

    if(use_bisection){
        MATRIX<T,TV::m> dKdU=h.dKdN*mK.dndU+Outer_Product(h.dKdphi,mK.N);
        MATRIX<T,TV::m> idKdK=((T)1-dKdU*mK.dUdX).Inverse();
        dYdK=dYdK*idKdK;
        dKdLi=dKdU*mK.dUdL;
        dKdAi=h.dKdN*mK.dndA+dKdU*mK.dUdA;}

    c.Y=Y;
    MATRIX<T,TV::m> dYdZ=dYdK*(h.dKdZ+h.dKdW*dWdZ);
    c.dYdLs=dYdZ*c.dZdLs;
    c.dYdAs=dYdZ*c.dZdAs;
    c.dYdLi=dYdK*(h.dKdX*dXdLi+h.dKdW*dWdLi+dKdLi)+dYdLi;
    c.dYdAi=dYdK*(h.dKdX*dXdAi+h.dKdW*dWdAi+dKdAi)+dYdAi;
}
//#####################################################################
// Function Update_Attachments_And_Prune_Pairs
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Update_Attachments_And_Prune_Pairs()
{
    TIMER_SCOPE_FUNC;
    int k=0;
    for(int i=0;i<collision_pairs.m;i++){
        COLLISION_PAIR c=collision_pairs(i);
        if(c.active){
            const RIGID_BODY<TV>& rbi=rigid_body_collection.Rigid_Body(c.bi);
            c.X=rbi.Frame().Inverse_Times(c.Y);
            collision_pairs(k++)=c;}
        else hash.Delete({c.bs,c.v,c.bi});}
    collision_pairs.Resize(k);
}
//#####################################################################
// Function Add_Pair
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Add_Pair(int bs,int v,int bi)
{
    TIMER_SCOPE_FUNC;
    // TODO: Interpolate X^n and X^(n+1) to choose surface point.
    if(hash.Contains({bs,v,bi})) return;
    const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(bs),
        &rbi=rigid_body_collection.Rigid_Body(bi);
    if(rbs.Has_Infinite_Inertia() && rbi.Has_Infinite_Inertia()) return;
    TV X=rbs.Frame()*rbs.simplicial_object->particles.X(v);
    if(rbi.implicit_object->Extended_Phi(X)>0) return;
    TV W=rbi.implicit_object->Closest_Point_On_Boundary(X);
    COLLISION_PAIR c={bs,v,bi,rbi.Frame().Inverse_Times(W)};
    collision_pairs.Append(c);
    hash.Insert({bs,v,bi});
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
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Read(TYPED_ISTREAM input)
{
    ARRAY<TRIPLE<int,int,int> > keys;
    Read_Binary(input,collision_pairs,keys);
    hash.Set_All(keys);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Write(TYPED_OSTREAM output) const
{
    ARRAY<TRIPLE<int,int,int> > keys;
    hash.Get_Keys(keys);
    keys.Sort();
    Write_Binary(output,collision_pairs,keys);
}
namespace PhysBAM{
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<float,2> >;
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<double,2> >;
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<float,3> >;
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<double,3> >;
}
