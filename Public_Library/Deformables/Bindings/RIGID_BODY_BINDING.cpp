//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Sergey Levine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Bindings/RIGID_BODY_BINDING.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_BINDING<TV>::
RIGID_BODY_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input)
    :BINDING<TV>(particles_input),rigid_body_collection(0),rigid_body_particles_index(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_BINDING<TV>::
RIGID_BODY_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input,const int particle_index_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
    const int rigid_body_particles_index_input,const TV& object_space_position_input)
    :BINDING<TV>(particles_input,particle_index_input),rigid_body_collection(&rigid_body_collection_input),rigid_body_particles_index(rigid_body_particles_index_input),
    object_space_position(object_space_position_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY_BINDING<TV>::
~RIGID_BODY_BINDING()
{
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV> RIGID_BODY_BINDING<TV>* RIGID_BODY_BINDING<TV>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    return new RIGID_BODY_BINDING(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(particles));
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
// TODO: This may not handle kinematic/static bodies
template<class TV> void RIGID_BODY_BINDING<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(particles.Size()+rigid_body_particles_index,particle_index));}
//#####################################################################
// Function Embedded_Position
//#####################################################################
template<class TV> TV RIGID_BODY_BINDING<TV>::
Embedded_Position() const
{
    return Rigid_Body().World_Space_Point(object_space_position);
}
//#####################################################################
// Function Embedded_Position
//#####################################################################
template<class TV> TV RIGID_BODY_BINDING<TV>::
Embedded_Position(ARRAY_VIEW<const TV> X) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV> TV RIGID_BODY_BINDING<TV>::
Embedded_Velocity() const
{
    return Rigid_Body().Pointwise_Object_Velocity(Embedded_Position());
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV> TV RIGID_BODY_BINDING<TV>::
Embedded_Velocity(ARRAY_VIEW<const TV> V) const
{
    PHYSBAM_FATAL_ERROR("Rigid body velocities required");
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV> TV RIGID_BODY_BINDING<TV>::
Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const
{
    return RIGID_BODY_STATE<TV>(rigid_body_collection->Rigid_Body(rigid_body_particles_index).Frame(),twist(rigid_body_particles_index)).Pointwise_Object_Velocity(Embedded_Position());
}
//#####################################################################
// Function Embedded_Acceleration
//#####################################################################
template<class TV> TV RIGID_BODY_BINDING<TV>::
Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench_input) const// TODO: consider including centrifugal acceleration
{
    RIGID_BODY<TV>& rigid_body=Rigid_Body();
    if(rigid_body.Has_Infinite_Inertia()) return TV();
    const TWIST<TV>& wrench=wrench_input(rigid_body_particles_index);
    return wrench.linear/rigid_body.Mass()+TV::Cross_Product(rigid_body.World_Space_Inertia_Tensor_Inverse()*wrench.angular,rigid_body.World_Space_Vector(object_space_position));
}
//#####################################################################
// Function One_Over_Effective_Mass
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_BODY_BINDING<TV>::
One_Over_Effective_Mass(const TV& direction) const // assumes direction is normalized
{
    RIGID_BODY<TV>& rigid_body=Rigid_Body();
    if(rigid_body.Has_Infinite_Inertia()) return 0;
    TV object_space_direction=rigid_body.Object_Space_Vector(direction);
    return TV::Dot_Product(object_space_direction,rigid_body.Object_Space_Impulse_Factor(object_space_position)*object_space_direction);
}
//#####################################################################
// Function One_Over_Effective_Mass
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_BODY_BINDING<TV>::
One_Over_Effective_Mass() const // return a lower bound for effective mass over all directions
{
    RIGID_BODY<TV>& rigid_body=Rigid_Body();if(rigid_body.Has_Infinite_Inertia()) return 0;return rigid_body.Object_Space_Impulse_Factor(object_space_position).Fast_Eigenvalues().Max();
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Apply_Impulse(const TV& impulse)
{
    Rigid_Body().Apply_Impulse_To_Body(Embedded_Position(),impulse);
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V) const
{
    PHYSBAM_ASSERT(Rigid_Body().Has_Infinite_Inertia());
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > rigid_V) const
{
    RIGID_BODY<TV>& rigid_body=Rigid_Body();
    if(rigid_body.Has_Infinite_Inertia()) return;
    rigid_V(rigid_body.particle_index).linear+=impulse/rigid_body.Mass();
    rigid_V(rigid_body.particle_index).angular+=rigid_body.World_Space_Inertia_Tensor_Inverse_Times(TV::Cross_Product(Rigid_Body().World_Space_Vector(object_space_position),impulse));
}
//#####################################################################
// Function Apply_Push
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Apply_Push(const TV& impulse)
{
    Rigid_Body().Apply_Push_To_Body(Embedded_Position(),impulse);
}
//#####################################################################
// Function Apply_Displacement_To_Parents_Based_On_Embedding
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Apply_Velocity_Change_To_Parents_Based_On_Embedding
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const
{
    PHYSBAM_FATAL_ERROR("Wrench required");
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const
{
    if(Rigid_Body().Has_Infinite_Inertia()) return;
    wrench_full(rigid_body_particles_index).linear+=force;
    wrench_full(rigid_body_particles_index).angular+=TV::Cross_Product(Rigid_Body().World_Space_Vector(object_space_position),force);
}
//#####################################################################
// Function Distribute_Mass_To_Parents
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const
{
    /*if(particles.mass(particle_index)) PHYSBAM_NOT_IMPLEMENTED();*/
}
//#####################################################################
// Function Impulse_Factor
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> RIGID_BODY_BINDING<TV>::
Impulse_Factor() const
{
    return Rigid_Body().Impulse_Factor(Rigid_Body().World_Space_Point(object_space_position));
}
//#####################################################################
// Function Parents
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Parents(ARRAY<int>& parents) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Weights
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Weights(ARRAY<T>& weights) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Read_Helper
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Read_Helper(TYPED_ISTREAM& input)
{
    BINDING<TV>::Read_Helper(input);Read_Binary(input,rigid_body_particles_index);Read_Binary(input,object_space_position);
}
//#####################################################################
// Function Write_Helper
//#####################################################################
template<class TV> void RIGID_BODY_BINDING<TV>::
Write_Helper(TYPED_OSTREAM& output) const
{
    BINDING<TV>::Write_Helper(output);Write_Binary(output,rigid_body_particles_index);Write_Binary(output,object_space_position);
}
namespace PhysBAM{
template class RIGID_BODY_BINDING<VECTOR<double,1> >;
template class RIGID_BODY_BINDING<VECTOR<double,2> >;
template class RIGID_BODY_BINDING<VECTOR<double,3> >;
template class RIGID_BODY_BINDING<VECTOR<float,1> >;
template class RIGID_BODY_BINDING<VECTOR<float,2> >;
template class RIGID_BODY_BINDING<VECTOR<float,3> >;
}
