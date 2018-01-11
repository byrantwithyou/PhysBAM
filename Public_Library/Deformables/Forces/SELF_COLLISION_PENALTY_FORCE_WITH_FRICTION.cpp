//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Matrices/IDENTITY_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/ZERO_MATRIX.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Forces/SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION(DEFORMABLE_PARTICLES<TV>& particles_input,
    T stiffness_coefficient,T friction)
    :BASE(particles_input),stiffness_coefficient(stiffness_coefficient),
    friction(friction)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
~SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active){
            const T_SURFACE& s=*surfaces(c.s);
            TV j=stiffness_coefficient*(particles.X(c.p)-c.Y);
            F(c.p)-=j;
            for(int k=0;k<TV::m;k++)
                F(s.mesh.elements(c.e)(k))+=c.w(k)*j;}}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    get_candidates();

    for(int i=0;i<collision_pairs.m;i++)
        Relax_Attachment(i);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=0;i<collision_pairs.m;i++){
        const COLLISION_PAIR& c=collision_pairs(i);
        if(c.active){
            VECTOR<int,TV::m> e=surfaces(c.s)->mesh.elements(c.e);
            TV j=stiffness_coefficient*(particles.X(c.p)-c.Y);
            TV dw=c.dwdZ*V(c.p),dY;
            for(int k=0;k<TV::m;k++)
                dw+=c.dwdE(k)*V(e(k));
            for(int k=0;k<TV::m;k++)
                dY+=c.w(k)*V(e(k))+dw(k)*particles.X(e(k));
            TV dj=stiffness_coefficient*(V(c.p)-dY);
            F(c.p)-=dj;
            for(int k=0;k<TV::m;k++)
                F(e(k))+=c.w(k)*dj+dw(k)*j;}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
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
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Relax_Attachment(int cp)
{
    // COLLISION_PAIR& c=collision_pairs(cp);
    // SELF_COLLISION<TV>* io=ios(c.o);
    // TV Z=particles.X(c.p);
    // T phi=io->Extended_Phi(Z);
    // TV n=io->Extended_Normal(Z);
    // SYMMETRIC_MATRIX<T,TV::m> H=io->Hessian(Z);
    // auto pr=Relax_Attachment_Helper(Z,c.X,phi,n,friction);
    // c.Y=pr.Y;
    // c.active=pr.active;
    // c.dYdZ=pr.dYdZ+Outer_Product(pr.dYdphi,n)+pr.dYdn*H;
}
//#####################################################################
// Function Update_Attachments_And_Prune_Pairs
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Attachments_And_Prune_Pairs()
{
    int k=0;
    for(int i=0;i<collision_pairs.m;i++){
        COLLISION_PAIR c=collision_pairs(i);
        if(c.active){
            c.w0=c.w;
            c.e0=c.e;
            collision_pairs(k++)=c;}
        else hash.Delete({c.p,c.s});}
    collision_pairs.Resize(k);
}
//#####################################################################
// Function Add_Pair
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Pair(int p,int s)
{
    // TODO: Interpolate X^n and X^(n+1) to choose surface point.
    if(hash.Contains({p,s})) return;
    COLLISION_PAIR c={p,s};
    TV Z=surfaces(s)->Closest_Point_On_Boundary(particles.X(p),trial_distance,0,&c.e0);
    auto elem=surfaces(s)->Get_Element(c.e0);
    c.w0=elem.Barycentric_Coordinates(Z);
    collision_pairs.Append(c);
    hash.Insert({p,s});
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
CFL_Strain_Rate() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
namespace PhysBAM{
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,2> >;
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<float,3> >;
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,2> >;
template class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<VECTOR<double,3> >;
}
