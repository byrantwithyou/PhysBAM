//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Matrices/MATRIX.h>
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
            F(c.p)-=stiffness_coefficient*V(c.p);}
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
    if(phi>=0){c.Y=Z;c.active=false;return;}
    TV n=io->Extended_Normal(Z);
    
    TV F=c.X-Z;
    T fn=F.Dot(n);
    if(fn<0){c.Y=Z;c.active=false;return;}
    T ft=(F-n*fn).Magnitude();
    c.active=true;
    if(ft<=friction*fn){c.Y=c.X;return;}

    // Project to friction cone.
    TV W=io->Closest_Point_On_Boundary(Z),H=W-c.X;
    T mu2=(1+sqr(friction)),Hn=H.Dot(n),Fn=F.Dot(n);
    QUADRATIC<T> q(mu2*sqr(Hn)-H.Dot(H),2*(mu2*Hn*Fn-H.Dot(F)),mu2*sqr(Fn)-F.Dot(F));
    q.Compute_Roots_In_Interval(0,1);
    PHYSBAM_ASSERT(q.roots==1);
    c.Y=io->Closest_Point_On_Boundary(q.root1*H+c.X);
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
}
