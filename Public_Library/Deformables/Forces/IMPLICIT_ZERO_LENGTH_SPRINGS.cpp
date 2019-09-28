//#####################################################################
// Copyright 2006-2007, Craig Schroeder, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_ZERO_LENGTH_SPRINGS
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/SPARSE_UNION_FIND.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/Robust_Arithmetic.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Deformables/Forces/IMPLICIT_ZERO_LENGTH_SPRINGS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using ::std::sqrt;
using namespace PhysBAM;
template<class TV> IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
IMPLICIT_ZERO_LENGTH_SPRINGS(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh_input)
    :DEFORMABLES_FORCES<TV>(particles),segment_mesh(segment_mesh_input)
{
    Set_Stiffness(0);Set_Damping(0);
}
template<class TV> IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
~IMPLICIT_ZERO_LENGTH_SPRINGS()
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    segment_mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    Update_Force_Elements(force_segments,segment_mesh.elements,particle_is_simulated);
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient) // assumes mass is already defined
{
    constant_stiffness=0;
    stiffness.Resize(segment_mesh.elements.m,no_init);
    for(int i=0;i<segment_mesh.elements.m;i++){
        int end1,end2;segment_mesh.elements(i).Get(end1,end2);
        T reduced_mass=Pseudo_Inverse(particles.one_over_effective_mass(end1)+particles.one_over_effective_mass(end2));
        stiffness(i)=scaling_coefficient*reduced_mass;}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;
    damping.Resize(segment_mesh.elements.m,no_init);
    for(int i=0;i<segment_mesh.elements.m;i++){
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(segment_mesh.elements(i)(0))+particles.one_over_effective_mass(segment_mesh.elements(i)(1)));
        T ym;if(!stiffness.m) ym=constant_stiffness;else ym=stiffness(i);
        damping(i)=overdamping_fraction*2*sqrt(ym*harmonic_mass);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;
    damping.Resize(segment_mesh.elements.m,no_init);
    for(int i=0;i<segment_mesh.elements.m;i++){
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(segment_mesh.elements(i)(0))+particles.one_over_effective_mass(segment_mesh.elements(i)(1)));
        T ym;if(!stiffness.m) ym=constant_stiffness;else ym=stiffness(i);
        damping(i)=overdamping_fraction(i)*2*sqrt(ym*harmonic_mass);}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    Add_Implicit_Velocity_Independent_Forces(particles.X,F,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    if(!damping.m) for(int s:force_segments){
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV force=constant_damping*(V(node1)-V(node2));
        F(node1)-=force;F(node2)+=force;}
    else for(int s:force_segments){
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV force=damping(s)*(V(node1)-V(node2));
        F(node1)-=force;F(node2)+=force;}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose) const
{
    if(!stiffness.m) for(int s:force_segments){
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV force=constant_stiffness*(V(node2)-V(node1));
        F(node1)+=force;F(node2)-=force;}
    else for(int s:force_segments){
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV force=stiffness(s)*(V(node2)-V(node1));
        F(node1)+=force;F(node2)-=force;}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(int s:force_segments){
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        T stiff=stiffness.m?stiffness(s):constant_stiffness;
        potential_energy+=(T).5*stiff*(particles.X(node2)-particles.X(node1)).Magnitude_Squared();}
    return potential_energy;
}
//#####################################################################
// Function Create_Edge_Zero_Length_Springs
//#####################################################################
template<class T,class TV> IMPLICIT_ZERO_LENGTH_SPRINGS<TV>* PhysBAM::
Create_Edge_Zero_Length_Springs(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const T stiffness,const T overdamping_fraction,const bool verbose)
{
    IMPLICIT_ZERO_LENGTH_SPRINGS<TV>* zls=new IMPLICIT_ZERO_LENGTH_SPRINGS<TV>(particles,segment_mesh);
    zls->Set_Stiffness(stiffness);
    zls->Set_Overdamping_Fraction(overdamping_fraction);
    return zls;
}
//#####################################################################
namespace PhysBAM{
template class IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<float,2> >;
template class IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<float,3> >;
template IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<float,3> >* Create_Edge_Zero_Length_Springs<float,VECTOR<float,3> >(DEFORMABLE_PARTICLES<VECTOR<float,3> >&,SEGMENT_MESH&,float,float,bool);
template class IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<double,2> >;
template class IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<double,3> >;
template IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<double,3> >* Create_Edge_Zero_Length_Springs<double,VECTOR<double,3> >(DEFORMABLE_PARTICLES<VECTOR<double,3> >&,SEGMENT_MESH&,double,double,bool);
}
