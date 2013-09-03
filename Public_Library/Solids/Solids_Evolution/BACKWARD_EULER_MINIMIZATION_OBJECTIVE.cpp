//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
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
    tmp1(static_cast<GENERALIZED_VELOCITY<TV>&>(*v1.Clone_Default())),minimization_system(minimization_system),collision_thickness(1e-8)
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
    tmp1=dv;
    Adjust_For_Collision(tmp1);
    Compute_Unconstrained(tmp1,h,g,e);
    if(g) Project_Gradient_And_Prune_Constraints(*g);
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Compute_Unconstrained(const KRYLOV_VECTOR_BASE<T>& Bdv,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const
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
// Function Adjust_For_Collision
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Adjust_For_Collision(KRYLOV_VECTOR_BASE<T>& Bdv) const
{
    GENERALIZED_VELOCITY<TV>& dv=debug_cast<GENERALIZED_VELOCITY<TV>&>(Bdv);
    minimization_system.colliding_particles.Remove_All();
    minimization_system.colliding_normals.Remove_All();
    if(!collision_objects.m) return;

    for(int p=0;p<dv.V.array.m;p++){
        TV dV=dv.V.array(p);
        TV V=v0.V.array(p)+dV;
        TV X=X0(p)+dt*V;
        T deepest_phi=collision_thickness;
        int deepest_index=-1;
        for(int j=0;j<collision_objects.m;j++){
            if(disabled_collision.Contains(PAIR<int,int>(j,p))) continue;
            IMPLICIT_OBJECT<TV>* io=collision_objects(j);
            T phi=io->Extended_Phi(X);
            if(phi<deepest_phi){
                deepest_phi=phi;
                deepest_index=j;}}
        if(deepest_index==-1) continue;
        IMPLICIT_OBJECT<TV>* io=collision_objects(deepest_index);
        for(int i=0;i<5 && abs(deepest_phi)>collision_thickness;i++){
            X-=deepest_phi*io->Extended_Normal(X);
            deepest_phi=io->Extended_Phi(X);}
        V=(X-X0(p))/dt;
        dV=V-v0.V.array(p);
        dv.V.array(p)=dV;
        minimization_system.colliding_particles.Append(p);
        minimization_system.colliding_normals.Append(io->Extended_Normal(X));}
}
//#####################################################################
// Function Project_Gradient_And_Prune_Constraints
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Project_Gradient_And_Prune_Constraints(KRYLOV_VECTOR_BASE<T>& Bg) const
{
    GENERALIZED_VELOCITY<TV>& g=debug_cast<GENERALIZED_VELOCITY<TV>&>(Bg);
    int sep=0;
    for(int i=minimization_system.colliding_particles.m-1;i>=0;i--){
        int p=minimization_system.colliding_particles(i);
        TV n=minimization_system.colliding_normals(i),&gv=g.V.array(p);
        T gv_n=gv.Dot(n);
        if(gv_n<0){
            sep++;
            minimization_system.colliding_particles.Remove_Index_Lazy(i);
            minimization_system.colliding_normals.Remove_Index_Lazy(i);}
        else gv-=gv_n*n;}
    if(minimization_system.example_forces_and_velocities)
        minimization_system.example_forces_and_velocities->Zero_Out_Enslaved_Velocity_Nodes(g.V.array,time,time);
}
//#####################################################################
// Function Make_Feasible
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Make_Feasible(KRYLOV_VECTOR_BASE<T>& dv) const
{
    Adjust_For_Collision(dv);
}
//#####################################################################
// Function Reset
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
//#####################################################################
// Function Disable_Current_Colliding_Pairs
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Disable_Current_Colliding_Pairs(T thickness)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    for(int j=0;j<collision_objects.m;j++){
        IMPLICIT_OBJECT<TV>* io=collision_objects(j);
        for(int p=0;p<particles.X.m;p++){
            T phi=io->Extended_Phi(particles.X(p));
            if(phi<=thickness) disabled_collision.Set(PAIR<int,int>(j,p));}}
}
//#####################################################################
// Function Test_Diff
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Test_Diff(const KRYLOV_VECTOR_BASE<T>& dv)
{
    KRYLOV_VECTOR_BASE<T> *t0=dv.Clone_Default();
    *t0=dv;
    Adjust_For_Collision(*t0);

    ARRAY<IMPLICIT_OBJECT<TV>*> collision_objects_copy;
    collision_objects_copy.Exchange(collision_objects);
    this->Test(*t0,minimization_system);
    collision_objects_copy.Exchange(collision_objects);

    T eps=(T)1e-6;
    static RANDOM_NUMBERS<T> random;
    KRYLOV_VECTOR_BASE<T> *ddv=dv.Clone_Default();
    KRYLOV_VECTOR_BASE<T> *g0=dv.Clone_Default();
    KRYLOV_VECTOR_BASE<T> *g1=dv.Clone_Default();
    KRYLOV_VECTOR_BASE<T> *a=dv.Clone_Default();
    KRYLOV_VECTOR_BASE<T> *b=dv.Clone_Default();

    T e0=0,e1=0;
    Compute(*t0,&minimization_system,g0,&e0);

    for(int i=0,n=ddv->Raw_Size();i<n;i++)
        ddv->Raw_Get(i)=random.Get_Uniform_Number(-eps,eps);
    minimization_system.Project(*ddv);
    b->Copy(1,*t0,*ddv);

    ARRAY<int> ai=minimization_system.colliding_particles;
    ARRAY<TV> an=minimization_system.colliding_normals;
    Adjust_For_Collision(*b);
    minimization_system.colliding_particles=ai;
    minimization_system.colliding_normals=an;
    ddv->Copy(-1,*t0,*b);

    minimization_system.Multiply(*ddv,*a);
    minimization_system.Project(*a);

    // new dv
    Compute(*b,&minimization_system,g1,&e1);
    minimization_system.Multiply(*ddv,*b);
    minimization_system.Project(*b);

    T test0a=(minimization_system.Inner_Product(*g1,*ddv)+minimization_system.Inner_Product(*g0,*ddv))/(2*eps);
    T test0b=(e1-e0)/eps;
    T test0=(test0b-test0a)/max(abs(test0a),(T)1e-20);

    a->Copy(1,*b,*a);
    a->Copy((T).5,*a);
    T test1a=sqrt(minimization_system.Inner_Product(*a,*a))/eps;

    b->Copy(-1,*g0,*g1);
    T test1b=sqrt(minimization_system.Inner_Product(*b,*b))/eps;

    a->Copy(-1,*b,*a);
    T test1=sqrt(minimization_system.Inner_Product(*a,*a))/(eps*max(abs(test1a),(T)1e-20));

    LOG::cout<<"RELATIVE ENERGY "<<-log10(abs(e1-e0)/abs(e0))<<std::endl;
    LOG::cout<<"energy diff test C "<<test0a<<"    "<<test0b<<"    "<<test0<<std::endl;
    LOG::cout<<"force diff test C "<<test1a<<"    "<<test1b<<"    "<<test1<<std::endl;

    delete ddv;
    delete g0;
    delete g1;
    delete a;
    delete b;
    delete t0;
}
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<float,1> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<float,2> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<float,3> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<double,1> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<double,2> >;
template class BACKWARD_EULER_MINIMIZATION_OBJECTIVE<VECTOR<double,3> >;
}
