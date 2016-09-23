//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_COLLISIONS
//#####################################################################
#include <Core/Matrices/FRAME.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>
#include <Solids/Collisions/BW_COLLISIONS.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/BW_BACKWARD_EULER_SYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BW_COLLISIONS<TV>::
BW_COLLISIONS(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input)
    :solid_body_collection(solid_body_collection_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BW_COLLISIONS<TV>::
~BW_COLLISIONS()
{}
//#####################################################################
// Function Detect_Cloth_Body_Contact
//#####################################################################
template<class TV> void BW_COLLISIONS<TV>::
Detect_Cloth_Body_Contact()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;

    // Detect cloth/body constraints
    for(int rb=0;rb<rigid_body_particles.Size();rb++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(rb);
        IMPLICIT_OBJECT<VECTOR<T,TV::m> >& object_space_implicit_object=*rigid_body.implicit_object->object_space_implicit_object;
        FRAME<TV> frame=rigid_body.Frame().Inverse();
        for(int p=0;p<particles.Size();p++)
            if(object_space_implicit_object.Lazy_Inside(frame*particles.X(p)))
                cloth_body_constraints.Append_Unique(PAIR<int,int>(p,rb));}
}
//#####################################################################
// Function Remove_Separating_Cloth_Body_Contacts
//#####################################################################
template<class TV> void BW_COLLISIONS<TV>::
Remove_Separating_Cloth_Body_Contacts(BW_BACKWARD_EULER_SYSTEM<TV>& system,KRYLOV_VECTOR_BASE<T>& R,KRYLOV_VECTOR_BASE<T>& B,KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& Q)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    R=B;system.Multiply(V,Q);R-=Q;
    VECTOR_T& actual_R=debug_cast<VECTOR_T&>(R);
    for(int i=cloth_body_constraints.m-1;i>=0;i--){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(cloth_body_constraints(i).y);
        int p=cloth_body_constraints(i).x;
        if(TV::Dot_Product(actual_R.V(p),rigid_body.Implicit_Geometry_Normal(particles.X(p)))<=0) // TODO make sure this is the right direction
            cloth_body_constraints.Remove_Index_Lazy(i);}
}
//#####################################################################
namespace PhysBAM{
template class BW_COLLISIONS<VECTOR<float,1> >;
template class BW_COLLISIONS<VECTOR<float,2> >;
template class BW_COLLISIONS<VECTOR<float,3> >;
template class BW_COLLISIONS<VECTOR<double,1> >;
template class BW_COLLISIONS<VECTOR<double,2> >;
template class BW_COLLISIONS<VECTOR<double,3> >;
}
