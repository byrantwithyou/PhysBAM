//#####################################################################
// Copyright 2006-2008, Ronald Fedkiw, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTILINEAR_SPRINGS
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/Robust_Arithmetic.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Deformables/Forces/MULTILINEAR_SPRINGS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Set_Spring_Phases
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Spring_Phases(const ARRAY<VECTOR<T,2> >& compression_intervals_input,const ARRAY<VECTOR<T,2> >& stretching_intervals_input)
{
    intervals.Resize(RANGE<VECTOR<int,1> >(VECTOR<int,1>(-compression_intervals_input.m),VECTOR<int,1>(stretching_intervals_input.m)),no_init);
    youngs_modulus_scaling.Resize(intervals.domain,no_init);
    correction_force_over_youngs_modulus.Resize(intervals.domain,no_init);
    intervals(0)=0;youngs_modulus_scaling(0)=1;correction_force_over_youngs_modulus(0)=0;
    for(int i=0;i<stretching_intervals_input.m;i++){ // stretching phases
        intervals(i)=stretching_intervals_input(i)[0];youngs_modulus_scaling(i)=stretching_intervals_input(i)[1];
        correction_force_over_youngs_modulus(i)=correction_force_over_youngs_modulus(i-1)+intervals(i)*(youngs_modulus_scaling(i-1)-youngs_modulus_scaling(i));}
    for(int i=0;i<compression_intervals_input.m;i++){ // compression phases
        intervals(-i)=compression_intervals_input(i)[0];youngs_modulus_scaling(-i)=compression_intervals_input(i)[1];
        correction_force_over_youngs_modulus(-i)=correction_force_over_youngs_modulus(1-i)+intervals(-i)*(youngs_modulus_scaling(1-i)-youngs_modulus_scaling(-i));}

    if(constant_youngs_modulus){
        constant_base_youngs_modulus=constant_youngs_modulus;constant_youngs_modulus=0;
        youngs_modulus.Resize(segment_mesh.elements.m);}
    else base_youngs_modulus=youngs_modulus;
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    states.Resize(segment_mesh.elements.m,no_init);
    current_lengths.Resize(segment_mesh.elements.m,no_init);
    correction_force.Resize(segment_mesh.elements.m,no_init);
    spring_count.Resize(intervals.domain,no_init);
    ARRAY_VIEW<const TV> X(particles.X);
    Invalidate_CFL();
    spring_count.Fill(0);
    for(int s:force_segments){
        typename BASE::STATE& state=states(s);
        const VECTOR<int,2>& nodes=segment_mesh.elements(s);
        state.nodes=nodes;
        state.direction=X(nodes[1])-X(nodes[0]);
        current_lengths(s)=state.direction.Normalize();
        T relative_deformation=(current_lengths(s)-visual_restlength(s))/restlength(s);
        int index=Find_Interval(relative_deformation);
        spring_count(index)++;
        T ym;if(constant_base_youngs_modulus) ym=constant_base_youngs_modulus;else ym=base_youngs_modulus(s);
        youngs_modulus(s)=youngs_modulus_scaling(index)*ym;
        correction_force(s)=correction_force_over_youngs_modulus(index)*youngs_modulus(s);
        damping(s)=springs_damping(index)(s);
        state.coefficient=damping(s)/restlength(s);}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    LINEAR_SPRINGS<TV>::Add_Velocity_Independent_Forces(F,time);
    for(int s:force_segments){
        const VECTOR<int,2>& nodes=segment_mesh.elements(s);
        TV force=correction_force(s)*states(s).direction;
        F(nodes[0])+=force;F(nodes[1])-=force;}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    springs_damping.Resize(intervals.domain,no_init);
    for(int i=springs_damping.domain.min_corner.x;i<springs_damping.domain.max_corner.x;i++){
        Set_All_Springs_To_Phase(i);
        LINEAR_SPRINGS<TV>::Set_Overdamping_Fraction(overdamping_fraction);
        springs_damping(i)=damping;}
}
//#####################################################################
// Function Set_Damping
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Damping(const T constant_damping_input)
{
    springs_damping.Resize(intervals.domain,no_init);
    for(int i=springs_damping.domain.min_corner.x;i<springs_damping.domain.max_corner.x;i++){
        springs_damping(i).Resize(segment_mesh.elements.m,no_init);
        springs_damping(i).Fill(constant_damping_input);}
}
//#####################################################################
// Function Set_Damping
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Damping(ARRAY_VIEW<const T> damping_input)
{
    springs_damping.Resize(intervals.domain,no_init);
    for(int i=springs_damping.domain.min_corner.x;i<springs_damping.domain.max_corner.x;i++) springs_damping(i)=damping_input;
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction) // 1 is critically damped
{
    springs_damping.Resize(intervals.domain,no_init);
    for(int i=springs_damping.domain.min_corner.x;i<springs_damping.domain.max_corner.x;i++){
        Set_All_Springs_To_Phase(i);
        LINEAR_SPRINGS<TV>::Set_Overdamping_Fraction(overdamping_fraction);
        springs_damping(i)=damping;}
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    for(int i=springs_damping.domain.min_corner.x;i<springs_damping.domain.max_corner.x;i++){
        Set_All_Springs_To_Phase(i);
        LINEAR_SPRINGS<TV>::Set_Overdamping_Fraction(overdamping_fraction);
        for(int k=0;k<damping.m;k++) springs_damping(i)(k)=max(springs_damping(i)(k),damping(k));}
}
//#####################################################################
// Function Set_All_Springs_To_Phase
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_All_Springs_To_Phase(const int phase_index) // 1 is critically damped
{
    if(constant_base_youngs_modulus) for(int s=0;s<segment_mesh.elements.m;s++) youngs_modulus(s)=youngs_modulus_scaling(phase_index)*constant_base_youngs_modulus;
    else for(int s=0;s<segment_mesh.elements.m;s++) youngs_modulus(s)=youngs_modulus_scaling(phase_index)*base_youngs_modulus(s);
}
//#####################################################################
namespace PhysBAM{
template class MULTILINEAR_SPRINGS<VECTOR<float,2> >;
template class MULTILINEAR_SPRINGS<VECTOR<float,3> >;
template class MULTILINEAR_SPRINGS<VECTOR<double,2> >;
template class MULTILINEAR_SPRINGS<VECTOR<double,3> >;
}
