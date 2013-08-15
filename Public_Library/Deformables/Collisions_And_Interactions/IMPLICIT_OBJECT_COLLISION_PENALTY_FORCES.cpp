//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{
using ::std::exp;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES(DEFORMABLE_PARTICLES<TV>& particles,IMPLICIT_OBJECT<TV>* implicit_object,T stiffness,T separation_parameter,T length_scale)
    :DEFORMABLES_FORCES<TV>(particles),implicit_object(implicit_object),own_implicit_object(false),stiffness(stiffness),separation_parameter(separation_parameter),
    length_scale(length_scale),pe(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
~IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES()
{
    if(own_implicit_object) delete implicit_object;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Update_Position_Based_State_Particle
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Position_Based_State_Particle(int p)
{
    T phi=implicit_object->Extended_Phi(particles.X(p));
    if(phi>separation_parameter) return;
    TV n=implicit_object->Extended_Normal(particles.X(p));
    T x=(separation_parameter-phi)/length_scale,exp_xm1=expm1(x);
    penetrating_particles.Append(p);
    pe+=stiffness*(exp_xm1-x);
    T a=stiffness*exp_xm1/length_scale,b=stiffness*(exp_xm1+1)/sqr(length_scale);
    grad_pe.Append(-a*n);
    H_pe.Append(b*SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(n)-a*implicit_object->Hessian(particles.X(p)));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    pe=0;
    penetrating_particles.Remove_All();
    grad_pe.Remove_All();
    H_pe.Remove_All();
    if(colliding_particles.m)
        for(int i=0;i<colliding_particles.m;i++)
            Update_Position_Based_State_Particle(colliding_particles(i));
    else
        for(int i=0;i<particles.X.m;i++)
            Update_Position_Based_State_Particle(i);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    F.Subset(penetrating_particles)-=grad_pe;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    Add_Implicit_Velocity_Independent_Forces(dX,dF,time);
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(int i=0;i<penetrating_particles.m;i++){
        int p=penetrating_particles(i);
        F(p)-=H_pe(i)*V(p);}
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    for(int i=0;i<penetrating_particles.m;i++){
        SYMMETRIC_MATRIX<T,TV::m>& H=H_pe(i);
        DIAGONAL_MATRIX<T,TV::m> D;
        MATRIX<T,TV::m> V;
        H.Fast_Solve_Eigenproblem(D,V);
        H=SYMMETRIC_MATRIX<T,TV::m>::Conjugate(V,D.Positive_Definite_Part());}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>::
Potential_Energy(const T time) const
{
    return pe;
}
template class IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<float,1> >;
template class IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<float,2> >;
template class IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<float,3> >;
template class IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<double,1> >;
template class IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<double,2> >;
template class IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<VECTOR<double,3> >;
}
