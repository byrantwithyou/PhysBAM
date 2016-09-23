//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOFT_BINDINGS
//#####################################################################
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Bindings/LINEAR_BINDING.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> LINEAR_BINDING<TV,d>::
LINEAR_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input)
    :BINDING<TV>(particles_input)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> LINEAR_BINDING<TV,d>::
LINEAR_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input,const int particle_index_input,const VECTOR<int,d>& parents_input,const VECTOR<T,d>& weights_input)
    :BINDING<TV>(particles_input,particle_index_input),parents(parents_input),weights(weights_input)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> LINEAR_BINDING<TV,d>::
LINEAR_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input,const int particle_index_input,const VECTOR<int,d>& parents_input,const VECTOR<T,d-1>& weights_input)
    :BINDING<TV>(particles_input,particle_index_input),parents(parents_input),weights(weights_input)
{
    STATIC_ASSERT(d>2); // this would be confusing in the segment case because interpolation fractions refer to the other vertex
    weights[d-1]=1-weights_input.Sum();
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> LINEAR_BINDING<TV,d>* LINEAR_BINDING<TV,d>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    return new LINEAR_BINDING(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(particles));
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV,int d> int LINEAR_BINDING<TV,d>::
Name() const
{
    return Static_Name();
}
//#####################################################################
// Function Embedded_Position
//#####################################################################
template<class TV,int d> TV LINEAR_BINDING<TV,d>::
Embedded_Position() const
{
    TV result;
    for(int i=0;i<d;i++) result+=weights[i]*particles.X(parents[i]);
    return result;
}
//#####################################################################
// Function Embedded_Position
//#####################################################################
template<class TV,int d> TV LINEAR_BINDING<TV,d>::
Embedded_Position(ARRAY_VIEW<const TV> X) const
{
    TV result;
    for(int i=0;i<d;i++) result+=weights[i]*X(parents[i]);
    return result;
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV,int d> TV LINEAR_BINDING<TV,d>::
Embedded_Velocity() const 
{
    TV result;
    for(int i=0;i<d;i++) result+=weights[i]*particles.V(parents[i]);
    return result;
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV,int d> TV LINEAR_BINDING<TV,d>::
Embedded_Velocity(ARRAY_VIEW<const TV> V) const 
{
    TV result;
    for(int i=0;i<d;i++) result+=weights[i]*V(parents[i]);
    return result;
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV,int d> TV LINEAR_BINDING<TV,d>::
Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const 
{
    TV result;
    for(int i=0;i<d;i++) result+=weights[i]*V(parents[i]);
    return result;
}
//#####################################################################
// Function Embedded_Acceleration
//#####################################################################
template<class TV,int d> TV LINEAR_BINDING<TV,d>::
Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench) const
{
    TV a;
    for(int i=0;i<d;i++) a+=weights[i]*particles.one_over_mass(parents[i])*F(parents[i]);
    return a;
}
//#####################################################################
// Function One_Over_Effective_Mass
//#####################################################################
template<class TV,int d> typename TV::SCALAR LINEAR_BINDING<TV,d>::
One_Over_Effective_Mass() const
{
    T result=0;
    for(int i=0;i<d;i++) result+=sqr(weights[i])*particles.one_over_mass(parents[i]);
    return result;
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Apply_Impulse(const TV& impulse)
{
    LINEAR_BINDING::Apply_Impulse(impulse,particles.V);
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V) const
{
    for(int i=0;i<d;i++) V(parents[i])+=particles.one_over_mass(parents[i])*weights[i]*impulse;
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > rigid_V) const
{
    LINEAR_BINDING::Apply_Impulse(impulse,particles.V);
}
//#####################################################################
// Function Apply_Push
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Apply_Push(const TV& impulse)
{
    for(int i=0;i<d;i++) particles.X(parents[i])+=particles.one_over_mass(parents[i])*weights[i]*impulse;
}
//#####################################################################
// Function Apply_Displacement_To_Parents_Based_On_Embedding
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle)
{
    VECTOR<T,d> c=weights/weights.Magnitude_Squared();
    for(int i=0;i<d;i++) if(!skip_particle || !(*skip_particle)(parents[i])) particles.X(parents[i])+=c[i]*dX;
}
//#####################################################################
// Function Apply_Velocity_Change_To_Parents_Based_On_Embedding
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle)
{
    VECTOR<T,d> c=weights/weights.Magnitude_Squared();
    for(int i=0;i<d;i++) if(!skip_particle || !(*skip_particle)(parents[i])) particles.V(parents[i])+=c[i]*dV;
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const
{
    for(int i=0;i<d;i++) F_full(parents[i])+=weights[i]*force;
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const
{
    for(int i=0;i<d;i++) F_full(parents[i])+=weights[i]*force;
}
//#####################################################################
// Function Distribute_Mass_To_Parents
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const
{
    for(int i=0;i<d;i++) mass_full(parents[i])+=weights[i]*particles.mass(particle_index);
}
//#####################################################################
// Function Parents
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Parents(ARRAY<int>& p) const
{
    p.Append_Elements(parents);
}
//#####################################################################
// Function Weights
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Weights(ARRAY<T>& w) const
{
    w.Append_Elements(weights);
}
//#####################################################################
// Function Impulse_Factor
//#####################################################################
template<class TV,int d> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> LINEAR_BINDING<TV,d>::
Impulse_Factor() const
{
    return SYMMETRIC_MATRIX<T,TV::m>()+LINEAR_BINDING::One_Over_Effective_Mass();
}
//#####################################################################
// Function Read_Helper
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Read_Helper(TYPED_ISTREAM& input)
{
    BINDING<TV>::Read_Helper(input);Read_Binary(input,parents,weights);
}
//#####################################################################
// Function Write_Helper
//#####################################################################
template<class TV,int d> void LINEAR_BINDING<TV,d>::
Write_Helper(TYPED_OSTREAM& output) const
{
    BINDING<TV>::Write_Helper(output);Write_Binary(output,parents,weights);
}
namespace PhysBAM{
template LINEAR_BINDING<VECTOR<double,1>,1>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,1> >&,int,VECTOR<int,1> const&,VECTOR<double,1> const&);
template LINEAR_BINDING<VECTOR<double,1>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,1> >&);
template LINEAR_BINDING<VECTOR<double,1>,3>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,1> >&);
template LINEAR_BINDING<VECTOR<double,1>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,1> >&);
template LINEAR_BINDING<VECTOR<double,2>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,2> >&);
template LINEAR_BINDING<VECTOR<double,2>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,2> >&,int,VECTOR<int,2> const&,VECTOR<double,2> const&);
template LINEAR_BINDING<VECTOR<double,2>,3>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,2> >&);
template LINEAR_BINDING<VECTOR<double,2>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,2> >&);
template LINEAR_BINDING<VECTOR<double,3>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,3> >&);
template LINEAR_BINDING<VECTOR<double,3>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,3> >&,int,VECTOR<int,2> const&,VECTOR<double,2> const&);
template LINEAR_BINDING<VECTOR<double,3>,3>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,3> >&);
template LINEAR_BINDING<VECTOR<double,3>,3>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,3> >&,int,VECTOR<int,3> const&,VECTOR<double,3> const&);
template LINEAR_BINDING<VECTOR<double,3>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,3> >&);
template LINEAR_BINDING<VECTOR<double,3>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,3> >&,int,VECTOR<int,4> const&,VECTOR<double,3> const&);
template LINEAR_BINDING<VECTOR<double,3>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<double,3> >&,int,VECTOR<int,4> const&,VECTOR<double,4> const&);
template LINEAR_BINDING<VECTOR<float,1>,1>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,1> >&,int,VECTOR<int,1> const&,VECTOR<float,1> const&);
template LINEAR_BINDING<VECTOR<float,1>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,1> >&);
template LINEAR_BINDING<VECTOR<float,1>,3>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,1> >&);
template LINEAR_BINDING<VECTOR<float,1>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,1> >&);
template LINEAR_BINDING<VECTOR<float,2>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,2> >&);
template LINEAR_BINDING<VECTOR<float,2>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,2> >&,int,VECTOR<int,2> const&,VECTOR<float,2> const&);
template LINEAR_BINDING<VECTOR<float,2>,3>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,2> >&);
template LINEAR_BINDING<VECTOR<float,2>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,2> >&);
template LINEAR_BINDING<VECTOR<float,3>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,3> >&);
template LINEAR_BINDING<VECTOR<float,3>,2>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,3> >&,int,VECTOR<int,2> const&,VECTOR<float,2> const&);
template LINEAR_BINDING<VECTOR<float,3>,3>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,3> >&);
template LINEAR_BINDING<VECTOR<float,3>,3>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,3> >&,int,VECTOR<int,3> const&,VECTOR<float,3> const&);
template LINEAR_BINDING<VECTOR<float,3>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,3> >&);
template LINEAR_BINDING<VECTOR<float,3>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,3> >&,int,VECTOR<int,4> const&,VECTOR<float,3> const&);
template LINEAR_BINDING<VECTOR<float,3>,4>::LINEAR_BINDING(DEFORMABLE_PARTICLES<VECTOR<float,3> >&,int,VECTOR<int,4> const&,VECTOR<float,4> const&);
template LINEAR_BINDING<VECTOR<double,1>,2>* LINEAR_BINDING<VECTOR<double,1>,2>::Create(GEOMETRY_PARTICLES<VECTOR<double,1> >&);
template LINEAR_BINDING<VECTOR<double,1>,3>* LINEAR_BINDING<VECTOR<double,1>,3>::Create(GEOMETRY_PARTICLES<VECTOR<double,1> >&);
template LINEAR_BINDING<VECTOR<double,1>,4>* LINEAR_BINDING<VECTOR<double,1>,4>::Create(GEOMETRY_PARTICLES<VECTOR<double,1> >&);
template LINEAR_BINDING<VECTOR<double,2>,2>* LINEAR_BINDING<VECTOR<double,2>,2>::Create(GEOMETRY_PARTICLES<VECTOR<double,2> >&);
template LINEAR_BINDING<VECTOR<double,2>,3>* LINEAR_BINDING<VECTOR<double,2>,3>::Create(GEOMETRY_PARTICLES<VECTOR<double,2> >&);
template LINEAR_BINDING<VECTOR<double,2>,4>* LINEAR_BINDING<VECTOR<double,2>,4>::Create(GEOMETRY_PARTICLES<VECTOR<double,2> >&);
template LINEAR_BINDING<VECTOR<double,3>,2>* LINEAR_BINDING<VECTOR<double,3>,2>::Create(GEOMETRY_PARTICLES<VECTOR<double,3> >&);
template LINEAR_BINDING<VECTOR<double,3>,3>* LINEAR_BINDING<VECTOR<double,3>,3>::Create(GEOMETRY_PARTICLES<VECTOR<double,3> >&);
template LINEAR_BINDING<VECTOR<double,3>,4>* LINEAR_BINDING<VECTOR<double,3>,4>::Create(GEOMETRY_PARTICLES<VECTOR<double,3> >&);
template LINEAR_BINDING<VECTOR<float,1>,2>* LINEAR_BINDING<VECTOR<float,1>,2>::Create(GEOMETRY_PARTICLES<VECTOR<float,1> >&);
template LINEAR_BINDING<VECTOR<float,1>,3>* LINEAR_BINDING<VECTOR<float,1>,3>::Create(GEOMETRY_PARTICLES<VECTOR<float,1> >&);
template LINEAR_BINDING<VECTOR<float,1>,4>* LINEAR_BINDING<VECTOR<float,1>,4>::Create(GEOMETRY_PARTICLES<VECTOR<float,1> >&);
template LINEAR_BINDING<VECTOR<float,2>,2>* LINEAR_BINDING<VECTOR<float,2>,2>::Create(GEOMETRY_PARTICLES<VECTOR<float,2> >&);
template LINEAR_BINDING<VECTOR<float,2>,3>* LINEAR_BINDING<VECTOR<float,2>,3>::Create(GEOMETRY_PARTICLES<VECTOR<float,2> >&);
template LINEAR_BINDING<VECTOR<float,2>,4>* LINEAR_BINDING<VECTOR<float,2>,4>::Create(GEOMETRY_PARTICLES<VECTOR<float,2> >&);
template LINEAR_BINDING<VECTOR<float,3>,2>* LINEAR_BINDING<VECTOR<float,3>,2>::Create(GEOMETRY_PARTICLES<VECTOR<float,3> >&);
template LINEAR_BINDING<VECTOR<float,3>,3>* LINEAR_BINDING<VECTOR<float,3>,3>::Create(GEOMETRY_PARTICLES<VECTOR<float,3> >&);
template LINEAR_BINDING<VECTOR<float,3>,4>* LINEAR_BINDING<VECTOR<float,3>,4>::Create(GEOMETRY_PARTICLES<VECTOR<float,3> >&);
}
