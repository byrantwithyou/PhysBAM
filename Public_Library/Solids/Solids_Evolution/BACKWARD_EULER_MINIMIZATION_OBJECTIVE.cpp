//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Bindings/BINDING_LIST.h>
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
    tmp1(static_cast<GENERALIZED_VELOCITY<TV>&>(*v1.Clone_Default())),minimization_system(minimization_system),collision_thickness(1e-15),last_energy(FLT_MAX),collisions_in_solve(true)
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
    if(h)
        minimization_system.forced_collisions.Clean_Memory();
    Make_Feasible(tmp1);
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(tmp1.V.array);
    Compute_Unconstrained(tmp1,h,g,e);
    if(h)
        for(int i=0;i<minimization_system.collisions.m;i++)
            minimization_system.collisions(i).phi=0;
    if(g) Project_Gradient_And_Prune_Constraints(*g,h);
    if(h){
        for(int i=0;i<minimization_system.collisions.m;i++){
            const COLLISION& c = minimization_system.collisions(i);
            minimization_system.forced_collisions.Insert(c.p,c.object);}}
    if(g){
        GENERALIZED_VELOCITY<TV>& gg=debug_cast<GENERALIZED_VELOCITY<TV>&>(*g);
        solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(gg.V.array);
        solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(gg.V.array);}
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
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions(particles.X);

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

    T energy=minimization_system.Inner_Product(dv,tmp0)/2;
    energy+=Update_Position_Based_State_Early_Out(time,true,h?FLT_MAX:last_energy);

    if(h) last_energy=energy;

    if(e) *e=energy;

    if(g){
        tmp1*=0;
        solid_body_collection.Add_Velocity_Independent_Forces(tmp1.V.array,tmp1.rigid_V.array,time);
        g->Copy(-dt,tmp1,tmp0);}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> typename TV::SCALAR BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Update_Position_Based_State_Early_Out(T time,bool is_position_update,T energy_early_out) const
{
    T energy=0;
    for(int k=0;k<solid_body_collection.solids_forces.m;k++)
        if(!solids_forces_lazy.Contains(k)){
            solid_body_collection.solids_forces(k)->Update_Position_Based_State(time);
            energy+=solid_body_collection.solids_forces(k)->Potential_Energy(time);}
    for(int k=0;k<solid_body_collection.rigid_body_collection.rigids_forces.m;k++)
        if(!rigids_forces_lazy.Contains(k)){
            solid_body_collection.rigid_body_collection.rigids_forces(k)->Update_Position_Based_State(time);
            energy+=solid_body_collection.rigid_body_collection.rigids_forces(k)->Potential_Energy(time);}
    for(int k=0;k<solid_body_collection.deformable_body_collection.deformables_forces.m;k++)
        if(!deformables_forces_lazy.Contains(k)){
            solid_body_collection.deformable_body_collection.deformables_forces(k)->Update_Position_Based_State(time,is_position_update);
            energy+=solid_body_collection.deformable_body_collection.deformables_forces(k)->Potential_Energy(time);}

    if(energy>energy_early_out) return energy;

    for(HASHTABLE<int>::ITERATOR it(solids_forces_lazy);it.Valid();it.Next()){
        solid_body_collection.solids_forces(it.Key())->Update_Position_Based_State(time);
        energy+=solid_body_collection.solids_forces(it.Key())->Potential_Energy(time);
        if(energy>energy_early_out) return energy;}
    for(HASHTABLE<int>::ITERATOR it(rigids_forces_lazy);it.Valid();it.Next()){
        solid_body_collection.rigid_body_collection.rigids_forces(it.Key())->Update_Position_Based_State(time);
        energy+=solid_body_collection.rigid_body_collection.rigids_forces(it.Key())->Potential_Energy(time);
        if(energy>energy_early_out) return energy;}
    for(HASHTABLE<int>::ITERATOR it(deformables_forces_lazy);it.Valid();it.Next()){
        solid_body_collection.deformable_body_collection.deformables_forces(it.Key())->Update_Position_Based_State(time,is_position_update);
        energy+=solid_body_collection.deformable_body_collection.deformables_forces(it.Key())->Potential_Energy(time);
        if(energy>energy_early_out) return energy;}
    return energy;
}
//#####################################################################
// Function Adjust_For_Collision
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Adjust_For_Collision(KRYLOV_VECTOR_BASE<T>& Bdv) const
{
    GENERALIZED_VELOCITY<TV>& dv=debug_cast<GENERALIZED_VELOCITY<TV>&>(Bdv);
    minimization_system.collisions.Remove_All();
    if(!collision_objects.m) return;

    for(int p=0;p<dv.V.array.m;p++){
        TV dV=dv.V.array(p);
        TV V=v0.V.array(p)+dV;
        TV X=X0(p)+dt*V;
        T deepest_phi=collision_thickness;
        int deepest_index=-1;
        if(int* object_index=minimization_system.forced_collisions.Get_Pointer(p)){
            deepest_index=*object_index;
            deepest_phi=collision_objects(deepest_index)->Extended_Phi(X);}
        for(int j=0;j<collision_objects.m;j++){
            if(disabled_collision.Contains(PAIR<int,int>(j,p))) continue;
            IMPLICIT_OBJECT<TV>* io=collision_objects(j);
            T phi=io->Extended_Phi(X);
            if(phi<deepest_phi){
                deepest_phi=phi;
                deepest_index=j;}}
        if(deepest_index==-1) continue;
        IMPLICIT_OBJECT<TV>* io=collision_objects(deepest_index);
        COLLISION c={deepest_index,p,deepest_phi,0,io->Extended_Normal(X),TV(),io->Hessian(X)};
        minimization_system.collisions.Append(c);
        for(int i=0;i<50 && abs(deepest_phi)>collision_thickness;i++){
            X-=deepest_phi*io->Extended_Normal(X);
            deepest_phi=io->Extended_Phi(X);}
        V=(X-X0(p))/dt;
        dV=V-v0.V.array(p);
        dv.V.array(p)=dV;}
}
//#####################################################################
// Function Project_Gradient_And_Prune_Constraints
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Project_Gradient_And_Prune_Constraints(KRYLOV_VECTOR_BASE<T>& Bg,bool allow_sep) const
{
    GENERALIZED_VELOCITY<TV>& g=debug_cast<GENERALIZED_VELOCITY<TV>&>(Bg);
    for(int i=minimization_system.collisions.m-1;i>=0;i--){
        COLLISION& c=minimization_system.collisions(i);
        TV &gv=g.V.array(c.p);
        c.n_dE=gv.Dot(c.n);
        if(allow_sep && c.n_dE<0) minimization_system.collisions.Remove_Index_Lazy(i);
        else{
            c.H_dE=c.H*gv;
            gv-=c.n_dE*c.n+c.phi*c.H_dE;}}
    if(minimization_system.example_forces_and_velocities)
        minimization_system.example_forces_and_velocities->Zero_Out_Enslaved_Velocity_Nodes(g.V.array,time,time);
}
//#####################################################################
// Function Make_Feasible
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Make_Feasible(KRYLOV_VECTOR_BASE<T>& dv) const
{
    if(collisions_in_solve)
        Adjust_For_Collision(dv);
}
//#####################################################################
// Function Initial_Guess
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_OBJECTIVE<TV>::
Initial_Guess(KRYLOV_VECTOR_BASE<T>& Bdv) const
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;
    GENERALIZED_VELOCITY<TV>& dv=debug_cast<GENERALIZED_VELOCITY<TV>&>(Bdv);
    T e0=0,e1=0;
    dv*=0;
    Compute(dv,0,&tmp0,&e0);
    for(int p=0;p<particles.number;p++)
        dv.V.array(p)=particles.mass(p)?tmp0.V.array(p)/-particles.mass(p):TV();
    for(int p=0;p<rigid_body_particles.number;p++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(p);
        if(rigid_body.Has_Infinite_Inertia()) continue;
        dv.rigid_V.array(p)=-rigid_body.Inertia_Inverse_Times(tmp0.rigid_V.array(p));}
    Compute(dv,0,0,&e1);
    if(e1>e0) dv*=0;
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
    Make_Feasible(*t0);

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

    ARRAY<COLLISION> ac;
    ac.Exchange(minimization_system.collisions);
    Make_Feasible(*b);
    ac.Exchange(minimization_system.collisions);
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
