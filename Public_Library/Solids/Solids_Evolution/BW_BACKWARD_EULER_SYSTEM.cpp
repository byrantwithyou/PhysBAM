//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_BACKWARD_EULER_SYSTEM
//#####################################################################
#include <Core/Matrices/MATRIX_3X3.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Solids/Collisions/BW_COLLISIONS.h>
#include <Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/BW_BACKWARD_EULER_SYSTEM.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BW_BACKWARD_EULER_SYSTEM<TV>::
BW_BACKWARD_EULER_SYSTEM(SOLID_BODY_COLLECTION<TV>& solid_body_collection,BW_COLLISIONS<TV>& bw_collisions_input,EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities,
    const T dt_input,const T time_input)
    :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(false,true),
    solid_body_collection(solid_body_collection),
    bw_collisions(bw_collisions_input),
    example_forces_and_velocities(example_forces_and_velocities),
    dt(dt_input),time(time_input),
    GV_F(solid_body_collection.New_Generalized_Velocity()),
    GV_B(solid_body_collection.New_Generalized_Velocity())
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BW_BACKWARD_EULER_SYSTEM<TV>::
~BW_BACKWARD_EULER_SYSTEM()
{
    delete &GV_F;
    delete &GV_B;
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
    Project_Helper(V,true);
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Force(const VECTOR_T& V,VECTOR_T& F) const
{
    F.V.array.Subset(solid_body_collection.deformable_body_collection.simulated_particles).Fill(TV());
    F.rigid_V.array.Subset(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles).Fill(TWIST<TV>());
    solid_body_collection.Add_Implicit_Velocity_Independent_Forces(V,F,time);
    F*=dt;
    solid_body_collection.Add_Velocity_Dependent_Forces(V,F,time);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF,bool transpose) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);VECTOR_T& F=debug_cast<VECTOR_T&>(BF);
    Force(V,F);
    for(int i=0;i<V.V.Size();i++) F.V(i)=V.V(i)-dt*solid_body_collection.deformable_body_collection.particles.one_over_mass(i)*F.V(i);
/*    for(int i=0;i<V.rigid_V.Size();i++) F.rigid_V(i)=V.rigid_V(i)-dt*(projection_data.mass.world_space_rigid_mass_inverse(i)*F.rigid_V(i));*/
}
//#####################################################################
// Function Project_Helper
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Project_Helper(KRYLOV_VECTOR_BASE<T>& BV,const bool negate) const
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

    // User constrained nodes
    VECTOR_T& V=debug_cast<VECTOR_T&>(BV);
    example_forces_and_velocities.Zero_Out_Enslaved_Velocity_Nodes(V.V.array,time,time);

    // Cloth/body contacts
    for(int i=0;i<bw_collisions.cloth_body_constraints.m;i++){
        int particle_index=bw_collisions.cloth_body_constraints(i).x;
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(bw_collisions.cloth_body_constraints(i).y);
        TV particle_location=particles.X(particle_index);
        TV normal=rigid_body.Implicit_Geometry_Normal(particle_location);
        TV relative_velocity=rigid_body.Pointwise_Object_Velocity(particle_location)-particles.V(particle_index);
        if(negate) V.V.array(particle_index)=TV::Dot_Product(relative_velocity,normal)*normal;
        else V.V.array(particle_index)=(MATRIX<T,TV::m>::Identity_Matrix()-MATRIX<T,TV::m>::Outer_Product(normal,normal))*V.V.array(particle_index);
    }
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
    Project_Helper(BV,false);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double BW_BACKWARD_EULER_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV0,const KRYLOV_VECTOR_BASE<T>& BV1) const
{
    const VECTOR_T& V0=debug_cast<const VECTOR_T&>(BV0),&V1=debug_cast<const VECTOR_T&>(BV1);
    double inner_product=0;
    for(int i=0;i<V0.V.Size();i++) inner_product+=V0.V(i).Dot(V1.V(i))*solid_body_collection.deformable_body_collection.particles.mass(i);
//+ARRAYS_COMPUTATIONS::Inner_Product_Double_Precision(projection_data.mass.world_space_rigid_mass,V0.rigid_V,V1.rigid_V);
    //TODO this should not be mass scaled
    return inner_product;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR BW_BACKWARD_EULER_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(BR);
    T convergence_norm_squared=sqr(R.V.Maximum_Magnitude());
/*    for(int p=0;p<R.rigid_V.Size();p++){
        const TWIST<TV>& twist=R.rigid_V(p);
        const RIGID_BODY_MASS<TV,true> &rigid_mass=projection_data.mass.world_space_rigid_mass(p),&rigid_mass_inverse=projection_data.mass.world_space_rigid_mass_inverse(p);
        convergence_norm_squared=max(convergence_norm_squared,
        twist.linear.Magnitude_Squared()+rigid_mass_inverse.mass*Dot_Product(twist.angular,rigid_mass.inertia_tensor*twist.angular));}*/
    T convergence_norm=sqrt(convergence_norm_squared);
    return convergence_norm;
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void BW_BACKWARD_EULER_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const
{
}
//#####################################################################
namespace PhysBAM{
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<float,1> >;
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<float,2> >;
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<float,3> >;
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<double,1> >;
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<double,2> >;
template class BW_BACKWARD_EULER_SYSTEM<VECTOR<double,3> >;
}
