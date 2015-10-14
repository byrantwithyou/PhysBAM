//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Deformables/Bindings/BINDING.h>
#include <Deformables/Bindings/LINEAR_BINDING_DYNAMIC.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LINEAR_BINDING_DYNAMIC<TV>::
LINEAR_BINDING_DYNAMIC(DEFORMABLE_PARTICLES<TV>& particles_input)
    :BINDING<TV>(particles_input)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LINEAR_BINDING_DYNAMIC<TV>::
LINEAR_BINDING_DYNAMIC(DEFORMABLE_PARTICLES<TV>& particles_input,const int particle_index_input,const int number_of_parents)
    :BINDING<TV>(particles_input,particle_index_input),parents(number_of_parents),weights(number_of_parents)
{}
//#####################################################################
// Function Create
//#####################################################################
template<class TV> LINEAR_BINDING_DYNAMIC<TV>* LINEAR_BINDING_DYNAMIC<TV>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    return new LINEAR_BINDING_DYNAMIC(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(particles));
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> int LINEAR_BINDING_DYNAMIC<TV>::
Name() const
{
    return Static_Name();
}
//#####################################################################
// Function Embedded_Position
//#####################################################################
template<class TV> TV LINEAR_BINDING_DYNAMIC<TV>::
Embedded_Position() const
{
    TV result;
    for(int i=0;i<parents.m;i++) result+=weights(i)*particles.X(parents(i));
    return result;
}
//#####################################################################
// Function Embedded_Position
//#####################################################################
template<class TV> TV LINEAR_BINDING_DYNAMIC<TV>::
Embedded_Position(ARRAY_VIEW<const TV> X) const
{
    TV result;
    for(int i=0;i<parents.m;i++) result+=weights(i)*X(parents(i));
    return result;
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV> TV LINEAR_BINDING_DYNAMIC<TV>::
Embedded_Velocity() const 
{
    TV result;
    for(int i=0;i<parents.m;i++) result+=weights(i)*particles.V(parents(i));
    return result;
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV> TV LINEAR_BINDING_DYNAMIC<TV>::
Embedded_Velocity(ARRAY_VIEW<const TV> V) const 
{
    TV result;
    for(int i=0;i<parents.m;i++) result+=weights(i)*V(parents(i));
    return result;
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV> TV LINEAR_BINDING_DYNAMIC<TV>::
Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const 
{
    TV result;
    for(int i=0;i<parents.m;i++) result+=weights(i)*V(parents(i));
    return result;
}
//#####################################################################
// Function Embedded_Acceleration
//#####################################################################
template<class TV> TV LINEAR_BINDING_DYNAMIC<TV>::
Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench) const
{
    TV a;
    for(int i=0;i<parents.m;i++) a+=weights(i)*particles.one_over_mass(parents(i))*F(parents(i));
    return a;
}
//#####################################################################
// Function One_Over_Effective_Mass
//#####################################################################
template<class TV> auto LINEAR_BINDING_DYNAMIC<TV>::
One_Over_Effective_Mass() const -> T
{
    T result=0;
    for(int i=0;i<parents.m;i++) result+=sqr(weights(i))*particles.one_over_mass(parents(i));
    return result;
}
//#####################################################################
// Function Impulse_Factor
//#####################################################################
template<class TV> auto LINEAR_BINDING_DYNAMIC<TV>::
Impulse_Factor() const -> SYMMETRIC_MATRIX<T,TV::m>
{
    return SYMMETRIC_MATRIX<T,TV::m>()+LINEAR_BINDING_DYNAMIC::One_Over_Effective_Mass();
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Apply_Impulse(const TV& impulse)
{
    LINEAR_BINDING_DYNAMIC::Apply_Impulse(impulse,particles.V);
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V) const
{
    for(int i=0;i<parents.m;i++) V(parents(i))+=particles.one_over_mass(parents(i))*weights(i)*impulse;
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > rigid_V) const
{
    LINEAR_BINDING_DYNAMIC::Apply_Impulse(impulse,particles.V);
}
//#####################################################################
// Function Apply_Push
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Apply_Push(const TV& impulse)
{
    for(int i=0;i<parents.m;i++) particles.X(parents(i))+=particles.one_over_mass(parents(i))*weights(i)*impulse;
}
//#####################################################################
// Function Apply_Displacement_To_Parents_Based_On_Embedding
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle)
{
    T one_over_weights_squared=1/weights.Magnitude_Squared();
    for(int i=0;i<parents.m;i++) if(!skip_particle || !(*skip_particle)(parents(i))) particles.X(parents(i))+=one_over_weights_squared*weights(i)*dX;
}
//#####################################################################
// Function Apply_Velocity_Change_To_Parents_Based_On_Embedding
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle)
{
    T one_over_weights_squared=1/weights.Magnitude_Squared();
    for(int i=0;i<parents.m;i++) if(!skip_particle || !(*skip_particle)(parents(i))) particles.V(parents(i))+=one_over_weights_squared*weights(i)*dV;
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const
{
    for(int i=0;i<parents.m;i++) F_full(parents(i))+=weights(i)*force;
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const
{
    for(int i=0;i<parents.m;i++) F_full(parents(i))+=weights(i)*force;
}
//#####################################################################
// Function Distribute_Mass_To_Parents
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const
{
    for(int i=0;i<parents.m;i++) mass_full(parents(i))+=weights(i)*particles.mass(particle_index);
}
//#####################################################################
// Function Parents
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Parents(ARRAY<int>& p) const
{
    p.Append_Elements(parents);
}
//#####################################################################
// Function Weights
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Weights(ARRAY<T>& w) const
{
    w.Append_Elements(weights);
}
//#####################################################################
// Function Read_Helper
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Read_Helper(TYPED_ISTREAM& input)
{
    BINDING<TV>::Read_Helper(input);
    Read_Binary(input,parents,weights);
}
//#####################################################################
// Function Write_Helper
//#####################################################################
template<class TV> void LINEAR_BINDING_DYNAMIC<TV>::
Write_Helper(TYPED_OSTREAM& output) const
{
    BINDING<TV>::Write_Helper(output);
    Write_Binary(output,parents,weights);
}
template class LINEAR_BINDING_DYNAMIC<VECTOR<double,1> >;
template class LINEAR_BINDING_DYNAMIC<VECTOR<double,2> >;
template class LINEAR_BINDING_DYNAMIC<VECTOR<double,3> >;
template class LINEAR_BINDING_DYNAMIC<VECTOR<float,1> >;
template class LINEAR_BINDING_DYNAMIC<VECTOR<float,2> >;
template class LINEAR_BINDING_DYNAMIC<VECTOR<float,3> >;
}
