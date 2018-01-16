//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/IDENTITY_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/ZERO_MATRIX.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input,
    T stiffness_coefficient,T friction)
    :BASE(particles_input),stiffness_coefficient(stiffness_coefficient),
    friction(friction)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
~IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active)
            F(c.p)-=stiffness_coefficient*(particles.X(c.p)-c.Y);}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    get_candidates();

    for(int i=0;i<collision_pairs.m;i++)
        Relax_Attachment(i);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active)
            F(c.p)-=stiffness_coefficient*(V(c.p)-c.dYdZ*V(c.p));}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Potential_Energy(const T time) const
{
    T pe=0;
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active)
            pe+=(T).5*stiffness_coefficient*(particles.X(c.p)-c.Y).Magnitude_Squared();}
    return pe;
}
// Z = colliding point, X = original attachment, Y = computed attachment
// phi = phi(Z), n = normal(Z)
namespace PhysBAM{
template<class TV,class T> RELAX_ATTACHMENT_HELPER<TV>
Relax_Attachment_Helper(const TV& Z,const TV& X,T phi,const TV& n,T mu)
{
    RELAX_ATTACHMENT_HELPER<TV> h;

    if(phi>=0){h.Y=Z;h.active=false;h.dYdZ+=1;return h;}
    TV F=X-Z;
    T Fn=F.Dot(n);
    if(Fn<0){h.Y=Z;h.active=false;h.dYdZ+=1;return h;}
    T qc=(1+sqr(mu))*sqr(Fn)-F.Dot(F);
    TV dqcdX=2*(1+sqr(mu))*Fn*n-2*F;
    TV dqcdZ=-dqcdX;
    TV dqcdn=2*(1+sqr(mu))*Fn*F;

    if(qc>=0){h.Y=X;h.active=true;h.dYdX+=1;return h;}
    T qb=Fn*phi*sqr(mu);
    TV dqbdX=n*phi*sqr(mu);
    TV dqbdZ=-dqbdX;
    TV dqbdn=F*phi*sqr(mu);
    T dqbdphi=Fn*sqr(mu);
    T qa=sqr(phi*mu);
    T dqadphi=2*phi*sqr(mu);

    // a^2*qc - 2*a*(1-a)*qb + (1-a)^2*qa = 0
    // exactly true: qc < 0, qb <= 0, qa >= 0; thus a unique in [0,1].
    T D=sqrt(sqr(qb)-qa*qc);
    TV dDdX=(2*qb*dqbdX-qa*dqcdX)/(2*D);
    TV dDdZ=(2*qb*dqbdZ-qa*dqcdZ)/(2*D);
    TV dDdn=(2*qb*dqbdn-qa*dqcdn)/(2*D);
    T dDdphi=(2*qb*dqbdphi-qc*dqadphi)/(2*D);

    T a=qa/(D+qb+qa);
    TV dadX=-qa/sqr(D+qb+qa)*(dqbdX+dDdX);
    TV dadZ=-qa/sqr(D+qb+qa)*(dqbdZ+dDdZ);
    TV dadn=-qa/sqr(D+qb+qa)*(dqbdn+dDdn);
    T dadphi=dqadphi/(D+qb+qa)-qa/sqr(D+qb+qa)*(dqbdphi+dqadphi+dDdphi);

    h.Y=(1-a)*(Z-n*phi)+a*X;
    h.dYdX=Outer_Product(Z-n*phi,-dadX)+Outer_Product(X,dadX)+a;
    h.dYdZ=Outer_Product(Z-n*phi,-dadZ)+(1-a)+Outer_Product(X,dadZ);
    h.dYdn=Outer_Product(Z-n*phi,-dadn)+(1-a)*(-phi)+Outer_Product(X,dadn);
    h.dYdphi=(Z-n*phi)*(-dadphi)+(1-a)*(-n)+dadphi*X;
    h.active=true;

    return h;
}
// Y = relaxed attachment, W = projected attachment on surface = Y-phi*n
// phi = phi(Y), n = normal(Y), H = Hessian(Y)
template<class TV,class T>
void Project_Attachment_To_Surface(
    TV& Y,
    MATRIX<T,TV::m>& dYdZ,
    T phi,const TV& n,const SYMMETRIC_MATRIX<T,TV::m>& H)
{
    Y-=phi*n;
    TV dWdphi=-n;
    DIAGONAL_MATRIX<T,TV::m> dWdn(-phi*(TV()+1));
    MATRIX<T,TV::m> dWdZ=dYdZ+Outer_Product(dWdphi,n)*dYdZ+dWdn*H*dYdZ;
    dYdZ=dWdZ;
}
}
//#####################################################################
// Function Relax_Attachment
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Relax_Attachment(int cp)
{
    COLLISION_PAIR& c=collision_pairs(cp);
    IMPLICIT_OBJECT<TV>* io=ios(c.o);
    TV Z=particles.X(c.p);
    T phi=io->Extended_Phi(Z);
    TV n=io->Extended_Normal(Z);
    SYMMETRIC_MATRIX<T,TV::m> H=io->Hessian(Z);
    auto pr=Relax_Attachment_Helper(Z,c.X,phi,n,friction);
    c.dYdZ=pr.dYdZ+Outer_Product(pr.dYdphi,n)+pr.dYdn*H;
    if(pr.active)
        Project_Attachment_To_Surface(
            pr.Y,
            c.dYdZ,
            io->Extended_Phi(pr.Y),
            io->Extended_Normal(pr.Y),
            io->Hessian(pr.Y));
    c.Y=pr.Y;
    c.active=pr.active;
}
//#####################################################################
// Function Update_Attachments_And_Prune_Pairs
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Attachments_And_Prune_Pairs()
{
    int k=0;
    for(int i=0;i<collision_pairs.m;i++){
        COLLISION_PAIR c=collision_pairs(i);
        if(c.active){
            c.X=c.Y;
            collision_pairs(k++)=c;}
        else hash.Delete({c.p,c.o});}
    collision_pairs.Resize(k);
}
//#####################################################################
// Function Add_Pair
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Pair(int p,int o)
{
    // TODO: Interpolate X^n and X^(n+1) to choose surface point.
    if(hash.Contains({p,o})) return;
    TV X=particles.X(p);
    if(ios(o)->Extended_Phi(X)>0) return;
    TV W=ios(o)->Closest_Point_On_Boundary(X);
    COLLISION_PAIR c={p,o,W};
    collision_pairs.Append(c);
    hash.Insert({p,o});
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
CFL_Strain_Rate() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
namespace PhysBAM{
template class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,2> >;
template class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,3> >;
template class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,2> >;
template class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,3> >;
template RELAX_ATTACHMENT_HELPER<VECTOR<double,2> > Relax_Attachment_Helper<VECTOR<double,2>,double>(
    VECTOR<double,2> const&,VECTOR<double,2> const&,double,VECTOR<double,2> const&,double);
template RELAX_ATTACHMENT_HELPER<VECTOR<double,3> > Relax_Attachment_Helper<VECTOR<double,3>,double>(
    VECTOR<double,3> const&,VECTOR<double,3> const&,double,VECTOR<double,3> const&,double);
template RELAX_ATTACHMENT_HELPER<VECTOR<float,2> > Relax_Attachment_Helper<VECTOR<float,2>,float>(
    VECTOR<float,2> const&,VECTOR<float,2> const&,float,VECTOR<float,2> const&,float);
template RELAX_ATTACHMENT_HELPER<VECTOR<float,3> > Relax_Attachment_Helper<VECTOR<float,3>,float>(
    VECTOR<float,3> const&,VECTOR<float,3> const&,float,VECTOR<float,3> const&,float);
template void Project_Attachment_To_Surface<VECTOR<double,2>,double>(
    VECTOR<double,2>&,MATRIX<double,2>&,double,const VECTOR<double,2>&,const SYMMETRIC_MATRIX<double,2>&);
template void Project_Attachment_To_Surface<VECTOR<double,3>,double>(
    VECTOR<double,3>&,MATRIX<double,3>&,double,const VECTOR<double,3>&,const SYMMETRIC_MATRIX<double,3>&);
template void Project_Attachment_To_Surface<VECTOR<float,2>,float>(
    VECTOR<float,2>&,MATRIX<float,2>&,float,const VECTOR<float,2>&,const SYMMETRIC_MATRIX<float,2>&);
template void Project_Attachment_To_Surface<VECTOR<float,3>,float>(
    VECTOR<float,3>&,MATRIX<float,3>&,float,const VECTOR<float,3>&,const SYMMETRIC_MATRIX<float,3>&);
}
