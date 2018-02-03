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
IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input)
    :BASE(particles_input)
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
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active){
            if(transpose)
                F(c.p)+=(c.dYdZ-1).Transpose_Times(stiffness_coefficient*V(c.p));
            else F(c.p)+=(c.dYdZ-1)*(stiffness_coefficient*V(c.p));}}
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
// W = Z projected to surface, mu = friction
namespace PhysBAM{
template<class TV,class T> void
Relax_Attachment_Helper(RELAX_ATTACHMENT_HELPER<TV>& h,const TV& Z,const TV& X,const TV& W,T mu)
{
    TV ZW=Z-W,XW=X-W;
    T zw=ZW.Normalize();
    T xw=XW.Normalize();
    T q=mu*zw;
    if(xw<=q){h.K=X;h.dKdX+=1;h.dynamic=false;/*LOG::puts("stick");*/return;}
    h.K=W+q*XW;
    h.dKdX=q/xw*((T)1-Outer_Product(XW,XW));
    h.dKdZ=Outer_Product(XW,mu*ZW);
    h.dKdW=(T)1-h.dKdX-h.dKdZ;
    h.dynamic=true;
    //LOG::puts("dynamic");
}
template<class TV,class T> void
Relax_Attachment_Helper_Search(RELAX_ATTACHMENT_HELPER<TV>& h,const TV& Z,const TV& X,const TV& W,const IMPLICIT_OBJECT<TV>* io,T mu)
{
    T mu_bar=1/(sqr(mu)+1);
    auto criterion=[X,W,Z,io,mu_bar](T s)
        {
            TV K=W+s*(X-W);
            T phi=io->Extended_Phi(K);
            TV N=io->Extended_Normal(K);
            TV Y=K-phi*N;
            TV u=Y-Z;
            return sqr(u.Dot(N))-mu_bar*u.Magnitude_Squared();
        };

    if(criterion(1)>=0){h.K=X;h.dKdX+=1;h.dynamic=false;/*puts("stick");*/return;}
    if(criterion(0)<=0){h.K=W;h.dKdW+=1;h.dynamic=true;/*puts("frictionless");*/return;}
    T hi=1,lo=0;
    while(hi-lo>std::numeric_limits<T>::epsilon()*2){
        T md=lo+(hi-lo)/2;
        if(criterion(md)>=0) lo=md;
        else hi=md;}

    T s=lo;
    TV U=X-W;
    TV K=W+s*U;
    T phi=io->Extended_Phi(K);
    TV N=io->Extended_Normal(K);
    TV Y=K-phi*N;
    TV u=Y-Z;
    T a=u.Dot(N);
    TV r=a*N-mu_bar*u;
    TV v=phi*r-a*u;
    T j=1/(r.Dot(U));
    
    h.K=K;
    h.dKdphi=a*j*(1-mu_bar)*U;
    h.dKdN=Outer_Product(U,j*v);
    h.dKdZ=Outer_Product(U,j*r);
    h.dKdW=(1-s)-(1-s)*h.dKdZ;
    h.dKdX=s-s*h.dKdZ;
    
    //puts("dynamic");
    h.dynamic=true;
}

// X = interior point, W = X projected to surface
template<class TV,class T>
bool Project_Attachment_To_Surface(TV& W,const IMPLICIT_OBJECT<TV>* io,const TV& X,
    MATRIX<T,TV::m>& dWdX,bool exit_if_sep,TV* nn=0,SYMMETRIC_MATRIX<T,TV::m>* HH=0)
{
    T phi=io->Extended_Phi(X);
    if(exit_if_sep && phi>0) return false;
    TV n=io->Extended_Normal(X);
    SYMMETRIC_MATRIX<T,TV::m> H=io->Hessian(X);
    W=X-phi*n;
    dWdX=(T)1-Outer_Product(n)-phi*H;
    if(nn) *nn=n;
    if(HH) *HH=H;
    return phi<=0;
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
    TV Z=particles.X(c.p),W,V,N;
    SYMMETRIC_MATRIX<T,TV::m> H;
    MATRIX<T,TV::m> dWdZ,dVdY;
    c.active=Project_Attachment_To_Surface(W,io,Z,dWdZ,true);
    if(!c.active) return;
    RELAX_ATTACHMENT_HELPER<TV> h;
    if(use_bisection) Relax_Attachment_Helper_Search(h,Z,c.X,W,io,friction);
    else Relax_Attachment_Helper(h,Z,c.X,W,friction);
    Project_Attachment_To_Surface(V,io,h.K,dVdY,false,&N,&H);
    MATRIX<T,TV::m> dKdZ=h.dKdZ+h.dKdW*dWdZ;
    if(use_bisection) dKdZ=((T)1-h.dKdN*H-Outer_Product(h.dKdphi,N)).Inverse()*dKdZ;
    c.Y=V;
    c.dYdZ=dVdY*dKdZ;
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
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Read(TYPED_ISTREAM& input)
{
    ARRAY<PAIR<int,int> > keys;
    Read_Binary(input,collision_pairs,keys);
    hash.Set_All(keys);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>::
Write(TYPED_OSTREAM& output) const
{
    ARRAY<PAIR<int,int> > keys;
    hash.Get_Keys(keys);
    keys.Sort();
    Write_Binary(output,collision_pairs,keys);
}
namespace PhysBAM{
template class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,2> >;
template class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,3> >;
template class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,2> >;
template class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,3> >;
template void Relax_Attachment_Helper<VECTOR<double,2>,double>(
    RELAX_ATTACHMENT_HELPER<VECTOR<double,2> >&,VECTOR<double,2> const&,
    VECTOR<double,2> const&,VECTOR<double,2> const&,double);
template void Relax_Attachment_Helper<VECTOR<float,2>,float>(
    RELAX_ATTACHMENT_HELPER<VECTOR<float,2> >&,VECTOR<float,2> const&,
    VECTOR<float,2> const&,VECTOR<float,2> const&,float);
template void Relax_Attachment_Helper<VECTOR<double,3>,double>(
    RELAX_ATTACHMENT_HELPER<VECTOR<double,3> >&,VECTOR<double,3> const&,
    VECTOR<double,3> const&,VECTOR<double,3> const&,double);
template void Relax_Attachment_Helper<VECTOR<float,3>,float>(
    RELAX_ATTACHMENT_HELPER<VECTOR<float,3> >&,VECTOR<float,3> const&,
    VECTOR<float,3> const&,VECTOR<float,3> const&,float);
template void Relax_Attachment_Helper_Search<VECTOR<double,2>,double>(
    RELAX_ATTACHMENT_HELPER<VECTOR<double,2> >&,VECTOR<double,2> const&,
    VECTOR<double,2> const&,VECTOR<double,2> const&,
    IMPLICIT_OBJECT<VECTOR<double,2> > const*,double);
template void Relax_Attachment_Helper_Search<VECTOR<float,2>,float>(
    RELAX_ATTACHMENT_HELPER<VECTOR<float,2> >&,VECTOR<float,2> const&,
    VECTOR<float,2> const&,VECTOR<float,2> const&,
    IMPLICIT_OBJECT<VECTOR<float,2> > const*,float);
template void Relax_Attachment_Helper_Search<VECTOR<double,3>,double>(
    RELAX_ATTACHMENT_HELPER<VECTOR<double,3> >&,VECTOR<double,3> const&,
    VECTOR<double,3> const&,VECTOR<double,3> const&,
    IMPLICIT_OBJECT<VECTOR<double,3> > const*,double);
template void Relax_Attachment_Helper_Search<VECTOR<float,3>,float>(
    RELAX_ATTACHMENT_HELPER<VECTOR<float,3> >&,VECTOR<float,3> const&,
    VECTOR<float,3> const&,VECTOR<float,3> const&,
    IMPLICIT_OBJECT<VECTOR<float,3> > const*,float);
}
