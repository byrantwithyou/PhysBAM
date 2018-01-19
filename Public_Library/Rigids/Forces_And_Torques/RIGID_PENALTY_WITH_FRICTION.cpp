//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
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
    const ARRAY<MOVE_RIGID_BODY_DIFF<TV> >& move_rb_diff,T stiffness_coefficient,
    T friction)
    :BASE(rigid_body_collection_input),
    stiffness_coefficient(stiffness_coefficient),friction(friction),
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
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(c.bs),
            &rbi=rigid_body_collection.Rigid_Body(c.bi);
        if(c.active){
            TV j=stiffness_coefficient*(c.Z-c.Y);
            rigid_F(c.bs)-=rbs.Gather(TWIST<TV>(j,typename TV::SPIN()),c.Z);
            rigid_F(c.bi)+=rbi.Gather(TWIST<TV>(j,typename TV::SPIN()),c.Y);}}
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
    ARRAY_VIEW<TWIST<TV> > rigid_F,const T time,bool transpose) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(c.bs),
            &rbi=rigid_body_collection.Rigid_Body(c.bi);
        if(c.active){
            TV j=stiffness_coefficient*(c.Z-c.Y);
            if(transpose){
                TV dLs=rigid_V(c.bs).linear,dLi=rigid_V(c.bi).linear;
                auto dAs=rigid_V(c.bs).angular,dAi=rigid_V(c.bi).angular;
                TV dZ=c.dZdLs*dLs+c.dZdAs*dAs;
                TV dY=c.dYdLs*dLs+c.dYdAs*dAs+c.dYdLi*dLi+c.dYdAi*dAi;
                TV dj=stiffness_coefficient*(dZ-dY);
                rigid_F(c.bs)-=rbs.Gather(TWIST<TV>(dj,(dZ-dLs).Cross(j)),c.Z);
                rigid_F(c.bi)+=rbi.Gather(TWIST<TV>(dj,(dY-dLi).Cross(j)),c.Y);}
            else{
                TWIST<TV> qs=rbs.Scatter(rigid_V(c.bs),c.Z);
                TWIST<TV> qi=rbi.Scatter(rigid_V(c.bi),c.Y);
                TV dj=stiffness_coefficient*(qi.linear-qs.linear);
                TV cs=j.Cross(qs.angular),ci=j.Cross(qi.angular);
                TV dY=ci-dj,dZ=dj-cs;
                rigid_F(c.bi).linear+=c.dYdLi.Transpose_Times(dY)-ci;
                rigid_F(c.bi).angular+=c.dYdAi.Transpose_Times(dY);
                rigid_F(c.bs).linear+=cs+c.dYdLs.Transpose_Times(dY)+c.dZdLs.Transpose_Times(dZ);
                rigid_F(c.bs).angular+=c.dYdAs.Transpose_Times(dY)+c.dZdAs.Transpose_Times(dZ);}}}
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
            const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(c.bs);
            TV X=rbs.Frame()*rbs.simplicial_object->particles.X(c.v);
            pe+=(T).5*stiffness_coefficient*(X-c.Y).Magnitude_Squared();}}
    return pe;
}
namespace PhysBAM{
//#####################################################################
// Function Project_Attachment_To_Surface
//#####################################################################
template<class TV,class T> bool
Project_Attachment_To_Surface(TV& W,const MOVE_RIGID_BODY_DIFF<TV>& mr,
    const IMPLICIT_OBJECT<TV>* io,const TV& X,MATRIX<T,TV::m>& dWdX,
    MATRIX<T,TV::m>& dWdL,MATRIX<T,TV::m,TV::SPIN::m>& dWdA,bool exit_if_sep)
{
    MATRIX<T,TV::m> dUdX,dUdL,dndN;
    MATRIX<T,TV::m,TV::SPIN::m> dUdA,dndA;
    TV U=mr.Frame_Inverse_Times(X,dUdX,dUdL,dUdA);
    T phi=io->Extended_Phi(U);
    if(exit_if_sep && phi>0) return false;
    TV N=io->Extended_Normal(U);
    TV n=mr.Rotate(N,dndN,dndA);
    SYMMETRIC_MATRIX<T,TV::m> dNdU=io->Hessian(U);
    W=X-phi*n;
    MATRIX<T,TV::m> dWdU=-phi*dndN*dNdU-Outer_Product(n,N);
    dWdX=(T)1+dWdU*dUdX;
    dWdL=dWdU*dUdL;
    dWdA=dWdU*dUdA-phi*dndA;
    return phi<=0;
}
}
//#####################################################################
// Function Relax_Attachment
//#####################################################################
template<class TV> void RIGID_PENALTY_WITH_FRICTION<TV>::
Relax_Attachment(int cp)
{
    COLLISION_PAIR& c=collision_pairs(cp);
    const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(c.bs),
        &rbi=rigid_body_collection.Rigid_Body(c.bi);
    const IMPLICIT_OBJECT<TV>* io=rbi.implicit_object->object_space_implicit_object;
    const MOVE_RIGID_BODY_DIFF<TV>& mrs=move_rb_diff(c.bs);
    const MOVE_RIGID_BODY_DIFF<TV>& mri=move_rb_diff(c.bi);
    MATRIX<T,TV::m> dXdv,dXdLi,dZdv,dWdZ,dWdLi,dVdY,dVdLi;
    MATRIX<T,TV::m,TV::SPIN::m> dXdAi,dZdAs,dWdAi,dVdAi;

    TV Xs=rbs.simplicial_object->particles.X(c.v);
    c.Z=mrs.Frame_Times(Xs,dZdv,c.dZdLs,c.dZdAs);
    TV X=mri.Frame_Times(c.X,dXdv,dXdLi,dXdAi),W,V;

    c.active=Project_Attachment_To_Surface(W,mri,io,c.Z,dWdZ,dWdLi,dWdAi,true);
    if(!c.active) return;

    auto pr=Relax_Attachment_Helper(c.Z,X,W,friction);

    Project_Attachment_To_Surface(V,mri,io,pr.Y,dVdY,dVdLi,dVdAi,false);
    c.Y=V;
    MATRIX<T,TV::m> dVdZ=dVdY*(pr.dYdZ+pr.dYdW*dWdZ);
    c.dYdLs=dVdZ*c.dZdLs;
    c.dYdAs=dVdZ*c.dZdAs;
    c.dYdLi=dVdY*(pr.dYdX*dXdLi+pr.dYdW*dWdLi)+dVdLi;
    c.dYdAi=dVdY*(pr.dYdX*dXdAi+pr.dYdW*dWdAi)+dVdAi;
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
    // TODO: Interpolate X^n and X^(n+1) to choose surface point.
    if(hash.Contains({bs,v,bi})) return;
    const RIGID_BODY<TV>& rbs=rigid_body_collection.Rigid_Body(bs),
        &rbi=rigid_body_collection.Rigid_Body(bi);
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
namespace PhysBAM{
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<float,2> >;
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<double,2> >;
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<float,3> >;
template class RIGID_PENALTY_WITH_FRICTION<VECTOR<double,3> >;
}
