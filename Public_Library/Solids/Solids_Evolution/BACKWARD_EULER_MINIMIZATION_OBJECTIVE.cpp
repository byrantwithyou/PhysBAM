//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_OBJECTIVE.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_SYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
BACKWARD_EULER_MINIMIZATION_OBJECTIVE(SOLID_BODY_COLLECTION<TV>& solid_body_collection,BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>& minimization_system)
    :solid_body_collection(solid_body_collection),
    v1(solid_body_collection.deformable_body_collection.particles.V,solid_body_collection.rigid_body_collection.rigid_body_particles.twist,solid_body_collection),
    v0(static_cast<GENERALIZED_VELOCITY<TV>&>(*v1.Clone_Default())),tmp0(static_cast<GENERALIZED_VELOCITY<TV>&>(*v1.Clone_Default())),
    tmp1(static_cast<GENERALIZED_VELOCITY<TV>&>(*v1.Clone_Default())),minimization_system(minimization_system)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
~BACKWARD_EULER_MINIMIZATION_OBJECTIVE()
{
    delete &v0;
    delete &tmp0;
    delete &tmp1;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Compute(const KRYLOV_VECTOR_BASE<T>& Bdv,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const
{
    const GENERALIZED_VELOCITY<TV>& dv=debug_cast<const GENERALIZED_VELOCITY<TV>&>(Bdv);
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;

    for(int p=0;p<particles.number;p++){
        TV dV=dv.V.array(p);
        TV& V1=v1.V.array(p);
        V1=v0.V.array(p)+dV;
        particles.X(p)=X0(p)+dt*V1;
        tmp0.V.array(p)=particles.mass(p)*dV;}

    for(int p=0;p<rigid_body_particles.number;p++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(p);
        if(rigid_body.Has_Infinite_Inertia()) continue;
        const TWIST<TV>& d_twist=dv.rigid_V.array(p);
        TWIST<TV>& twist=v1.rigid_V.array(p);
        twist=v0.rigid_V.array(p)+d_twist;
        FRAME<TV>& frame=rigid_body_particles.frame(p);
        frame.t=frame0(p).t+dt*twist.linear;
        frame.r=ROTATION<TV>::From_Rotation_Vector(dt*twist.angular)*frame0(p).r;
        tmp0.rigid_V.array(p)=rigid_body.Inertia_Times(d_twist);}
    solid_body_collection.Update_Position_Based_State(time,true);

    if(e){
        T ke=0,pe=0;
        solid_body_collection.Compute_Energy(time,ke,pe);
        *e=minimization_system.Inner_Product(dv,tmp0)/2+pe;}

    if(g){
        tmp1*=0;
        solid_body_collection.Add_Velocity_Independent_Forces(tmp1.V.array,tmp1.rigid_V.array,time);
        g->Copy(-dt,tmp1,tmp0);}
}
//#####################################################################
// Function Resize_Vectors
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Reset()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;

    X0=particles.X;
    frame0=rigid_body_particles.frame;

    v1.V.array.Set(particles.V);
    v1.rigid_V.array.Set(rigid_body_particles.twist);
    v1.kinematic_and_static_rigid_V.array.Set(rigid_body_particles.twist);
    v0.Resize(v1);
    v0=v1;
    tmp0.Resize(v1);
    tmp1.Resize(v1);
}
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<float,1> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<float,2> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<float,3> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<double,1> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<double,2> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<double,3> >;
}
