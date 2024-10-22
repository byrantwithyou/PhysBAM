//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Log/FINE_TIMER.h>
#include <Core/Matrices/MATRIX.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/RELAX_ATTACHMENT_IMPLICIT.h>
#include <Rigids/Collisions/RELAX_ATTACHMENT_MESH.h>
#include <Rigids/Forces_And_Torques/MOVE_RIGID_BODY_DIFF.h>
#include <Rigids/Forces_And_Torques/RIGID_PENALTY_WITH_FRICTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/RIGID_DEFORMABLE_PENALTY_WITH_FRICTION.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
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
Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    TIMER_SCOPE_FUNC;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(c.b);
        if(c.active){
            TV j=stiffness_coefficient*(particles.X(c.p)-c.Y);
            F.V.array(c.p)-=j;
            if(!rb.Has_Infinite_Inertia())
                F.rigid_V.array(c.b)+=rb.Gather(TWIST<TV>(j,typename TV::SPIN()),c.Y);}}
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
Add_Implicit_Velocity_Independent_Forces(const GENERALIZED_VELOCITY<TV>& V,
    GENERALIZED_VELOCITY<TV>& F,const T time,bool transpose) const
{
    TIMER_SCOPE_FUNC;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(c.b);
        if(c.active){
            TV j=stiffness_coefficient*(particles.X(c.p)-c.Y);
            if(rb.Has_Infinite_Inertia()){
                if(transpose)
                    F.V.array(c.p)+=(c.dYdZ-1).Transpose_Times(stiffness_coefficient*V.V.array(c.p));
                else F.V.array(c.p)+=(c.dYdZ-1)*(stiffness_coefficient*V.V.array(c.p));}
            else{
                if(transpose){
                    TWIST<TV> tw=rb.Scatter(V.rigid_V.array(c.b),c.Y);
                    TV dj=tw.linear-V.V.array(c.p),cp=j.Cross(tw.angular);
                    TV dZ=stiffness_coefficient*dj,dY=cp-dZ;
                    F.rigid_V.array(c.b).linear+=c.dYdL.Transpose_Times(dY)-cp;
                    F.rigid_V.array(c.b).angular+=c.dYdA.Transpose_Times(dY);
                    F.V.array(c.p)+=dZ+c.dYdZ.Transpose_Times(dY);}
                else{
                    TV dZ=V.V.array(c.p),dL=V.rigid_V.array(c.b).linear;
                    auto dA=V.rigid_V.array(c.b).angular;
                    TV dY=c.dYdZ*dZ+c.dYdL*dL+c.dYdA*dA;
                    TV dj=stiffness_coefficient*(dZ-dY);
                    F.V.array(c.p)-=dj;
                    F.rigid_V.array(c.b)+=rb.Gather(TWIST<TV>(dj,(dY-dL).Cross(j)),c.Y);}}}}
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
    if(collision_pairs(cp).e0>=0) Relax_Attachment_Mesh(cp);
    else Relax_Attachment_Implicit(cp);
}
//#####################################################################
// Function Relax_Attachment_Implicit
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Relax_Attachment_Implicit(int cp)
{
    TIMER_SCOPE_FUNC;
    COLLISION_PAIR& c=collision_pairs(cp);
    const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(c.b);
    const IMPLICIT_OBJECT<TV>* io=rb.implicit_object->object_space_implicit_object;
    const MOVE_RIGID_BODY_DIFF<TV>& mr=move_rb_diff(c.b);
    MATRIX<T,TV::m> dXdv,dXdL,dZdv,dWdZ,dWdL,dYdK,dYdL,dKdL;
    MATRIX<T,TV::m,TV::SPIN::m> dXdA,dWdA,dYdA,dKdA;

    TV Z=particles.X(c.p);
    TV X=mr.Frame_Times(c.X,dXdv,dXdL,dXdA),W;

    MOVING_LEVEL_SET_HELPER<TV> mZ;
    c.active=mZ.Init(mr,io,Z,true);
    if(!c.active) return;
    Project_Attachment_To_Surface(W,mZ,Z,dWdZ,dWdL,dWdA);

    RELAX_ATTACHMENT_IMPLICIT<TV> h;
    if(use_bisection) h.Relax_Search(Z,X,W,rb.implicit_object,friction);
    else h.Relax(Z,X,W,friction);
    if(h.dynamic) num_dynamic++;
    else num_stick++;
    
    MOVING_LEVEL_SET_HELPER<TV> mK;
    mK.Init(mr,io,h.K,false);
    Project_Attachment_To_Surface(c.Y,mK,h.K,dYdK,dYdL,dYdA);

    if(use_bisection){
        MATRIX<T,TV::m> dKdU=h.dKdN*mK.dndU+Outer_Product(h.dKdphi,mK.N);
        MATRIX<T,TV::m> idKdK=((T)1-dKdU*mK.dUdX).Inverse();
        dYdK=dYdK*idKdK;
        dKdL=dKdU*mK.dUdL;
        dKdA=h.dKdN*mK.dndA+dKdU*mK.dUdA;}
    
    MATRIX<T,TV::m> dYdZ=dYdK*(h.dKdZ+h.dKdW*dWdZ);
    c.dYdZ=dYdZ;
    c.dYdL=dYdK*(h.dKdX*dXdL+h.dKdW*dWdL+dKdL)+dYdL;
    c.dYdA=dYdK*(h.dKdX*dXdA+h.dKdW*dWdA+dKdA)+dYdA;
}
//#####################################################################
// Function Relax_Attachment_Mesh
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Relax_Attachment_Mesh(int cp)
{
    TIMER_SCOPE_FUNC;
    COLLISION_PAIR& c=collision_pairs(cp);
    const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(c.b);
    const MOVE_RIGID_BODY_DIFF<TV>& mr=move_rb_diff(c.b);

    MATRIX<T,TV::m> dZidZ,dZidL,dYdYi,dYdL;
    MATRIX<T,TV::m,TV::SPIN::m> dZidA,dYdA;

    TV Z=particles.X(c.p);
    TV Zi=mr.Frame_Inverse_Times(Z,dZidZ,dZidL,dZidA);
    RELAX_ATTACHMENT_MESH<TV> ram;
    ram.Relax(c.e0,c.w0,Zi,*rb.simplicial_object,-1,friction);
    c.active=ram.active;
    if(!c.active) return;

    c.Y=mr.Frame_Times(ram.Y,dYdYi,dYdL,dYdA);
    c.e=ram.e;
    c.w=ram.w;

    MATRIX<T,TV::m> dYidL,dYidZ;
    MATRIX<T,TV::m,TV::SPIN::m> dYidA;

    for(int i=0;i<ram.diff_entry.m;i++){
        const auto& de=ram.diff_entry(i);
        dYidZ=de.dYdI(0)*dYidZ+de.dYdI(1)*dZidZ;
        dYidL=de.dYdI(0)*dYidL+de.dYdI(1)*dZidL;
        dYidA=de.dYdI(0)*dYidA+de.dYdI(1)*dZidA;}

    c.dYdZ=dYdYi*dYidZ;
    c.dYdL=dYdYi*dYidL+dYdL;
    c.dYdA=dYdYi*dYidA+dYdA;
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
            c.e0=c.e;
            c.w0=c.w;
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
    c.e0=c.e=-1;
    collision_pairs.Append(c);
    hash.Insert({p,b});
}
//#####################################################################
// Function Add_Pair
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Pair(int p,int b,int e,const TV& X0,const FRAME<TV>& f,T thickness)
{
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    TIMER_SCOPE_FUNC;
    // TODO: Interpolate X^n and X^(n+1) to choose surface point.
    if(hash.Contains({p,b})) return;
    TV X1=particles.X(p);
    const RIGID_BODY<TV>& rb=rigid_body_collection.Rigid_Body(b);
    const auto& ts=*rb.simplicial_object;
    TV Xb1=rb.Frame().Inverse_Times(X1);
    if(ts.Get_Element(e).Signed_Distance(Xb1)>0) return;
    VECTOR<TV,TV::m> E(ts.particles.X.Subset(ts.mesh.elements(e))),E0,E1;
    for(int i=0;i<TV::m;i++){
        E0(i)=f*E(i);
        E1(i)=rb.Frame()*E(i);}

    // Do the expensive check.
    T_FACE face(E0);
    T collision_time=0;
    TV normal;
    VECTOR<T,TV::m+1> weights;
    VECTOR<TV,TV::m> V_f(E1-E0);
    bool in=face.Point_Face_Collision(X0,X1-X0,V_f,1,
        thickness,collision_time,normal,weights,false);
    if(!in) return;

    TV w=weights.Remove_Index(0);
    if(w.Min()<0) return;
    TV W=E.Weighted_Sum(w);

    COLLISION_PAIR c={p,b,W};
    c.w0=w;
    c.e0=e;
    collision_pairs.Append(c);
    hash.Insert({p,b});
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces(const GENERALIZED_VELOCITY<TV>& V,
    GENERALIZED_VELOCITY<TV>& F,const T time) const
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
Add_Velocity_Dependent_Forces_First_Half(const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,GENERALIZED_VELOCITY<TV>& F,const T time) const
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
Read(TYPED_ISTREAM input)
{
    ARRAY<PAIR<int,int> > keys;
    Read_Binary(input,collision_pairs,keys);
    hash.Set_All(keys);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>::
Write(TYPED_OSTREAM output) const
{
    ARRAY<PAIR<int,int> > keys;
    hash.Get_Keys(keys);
    keys.Sort();
    Write_Binary(output,collision_pairs,keys);
}
namespace PhysBAM{
template class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<VECTOR<float,2> >;
template class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<VECTOR<double,2> >;
template class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<VECTOR<float,3> >;
template class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<VECTOR<double,3> >;
}
