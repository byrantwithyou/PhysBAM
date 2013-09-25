//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES(DEFORMABLE_PARTICLES<TV>& particles,TETRAHEDRALIZED_VOLUME<T>& collision_body,ARRAY<T>& undeformed_phi,
        T stiffness,T separation_parameter,T length_scale)
    :DEFORMABLES_FORCES<TV>(particles),collision_body(collision_body),undeformed_phi(undeformed_phi),stiffness(stiffness),separation_parameter(separation_parameter),
    length_scale(length_scale),pe(0)
{
    collision_body.Initialize_Hierarchy();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
~DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES()
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Update_Position_Based_State_Particle
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Position_Based_State_Particle(int p)
{
    TETRAHEDRON_MESH& tetrahedron_mesh=collision_body.mesh;
    ARRAY<int> particles_to_ignore;
    particles_to_ignore.Append(p);
    TV x=particles.X(p);
    ARRAY<int> intersection_list;
    collision_body.hierarchy->Intersection_List(x,intersection_list);
    for(int n=0;n<intersection_list.m;n++){
        int t=intersection_list(n);
        VECTOR<int,4> nodes=tetrahedron_mesh.elements(t);
        if(nodes.Contains(p)) continue;
        int i,j,k,l;
        nodes.Get(i,j,k,l);
        VECTOR<T,TV::m+1> w=TETRAHEDRON<T>::Barycentric_Coordinates(x,particles.X.Subset(nodes));
        if(w.Min()<-1e-12) continue;
        if(undeformed_phi.Subset(nodes).Weighted_Sum(w)>separation_parameter) continue;
        penetrating_particles.Append(nodes.Insert(p,0));
        auto X=From_Var<4,0>(particles.X(p)-particles.X(l));
        auto A=From_Var<4,1>(particles.X(i)-particles.X(l));
        auto B=From_Var<4,2>(particles.X(j)-particles.X(l));
        auto C=From_Var<4,3>(particles.X(k)-particles.X(l));
        auto B_cross_C=B.Cross(C);
        auto Det=A.Dot(B_cross_C);
        auto Y=X/Det;
        auto wA=Y.Dot(B_cross_C);
        auto C_cross_A=C.Cross(A);
        auto wB=Y.Dot(C_cross_A); 
        auto A_cross_B=A.Cross(B);
        auto wC=Y.Dot(A_cross_B); 
        auto wD=1-(wA+wB+wC);
        auto phi=undeformed_phi(i)*wA+undeformed_phi(j)*wB+undeformed_phi(k)*wC+undeformed_phi(l)*wD;
        auto s=(separation_parameter-phi)/length_scale;
        auto ee=stiffness*(exp(s)-s-1);
        pe+=ee.x;
        VECTOR<TV,5> de;
        for(int i=0;i<4;i++){
            TV t=ee.dx(i);
            de(i)=t;
            de(4)-=t;}
        grad_pe.Append(de);
        VECTOR<VECTOR<MATRIX<T,TV::m>,5>,5> he; 
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                MATRIX<T,TV::m> t=ee.ddx(i,j);
                he(i)(j)=t;
                he(i)(4)-=t;
                he(4)(4)+=t;}
            he(4)(i)=he(i)(4).Transposed();}
        H_pe.Append(he);
        break;}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    pe=0;
    penetrating_particles.Remove_All();
    grad_pe.Remove_All();
    H_pe.Remove_All();
    collision_body.hierarchy->Update_Boxes(separation_parameter);
    if(colliding_particles.m)
        for(int pp=0;pp<colliding_particles.m;pp++){
            int p=colliding_particles(pp);
            Update_Position_Based_State_Particle(p);}
    else
        for(int p=0;p<particles.X.m;p++){
            Update_Position_Based_State_Particle(p);}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(int pp=0;pp<penetrating_particles.m;pp++)
        F.Subset(penetrating_particles(pp))-=grad_pe(pp);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    Add_Implicit_Velocity_Independent_Forces(dX,dF,time);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(int pp=0;pp<penetrating_particles.m;pp++){
        const VECTOR<int,TV::m+2>& n=penetrating_particles(pp);
        for(int i=0;i<n.m;i++){
            int p=n(i);
            for(int j=0;j<n.m;j++)
                F(p)-=H_pe(pp)(i)(j)*V(n(j));}}
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Potential_Energy(const T time) const
{
    return pe;
}
template class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<float,3> >;
template class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<double,3> >;
}
