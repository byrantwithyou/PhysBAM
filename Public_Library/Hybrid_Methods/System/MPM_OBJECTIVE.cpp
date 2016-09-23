//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_SYSTEM.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
#include <Hybrid_Methods/System/MPM_OBJECTIVE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_OBJECTIVE<TV>::
MPM_OBJECTIVE(MPM_EXAMPLE<TV>& example)
    :system(*new MPM_KRYLOV_SYSTEM<TV>(example)),
    v0(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    v1(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    tmp0(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    tmp1(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    tmp2(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices)),
    collision_thickness(1e-15)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_OBJECTIVE<TV>::
~MPM_OBJECTIVE()
{
    delete &system;
    delete &v0;
    delete &v1;
    delete &tmp0;
    delete &tmp1;
    delete &tmp2;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Compute(const KRYLOV_VECTOR_BASE<T>& Bdv,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const
{
    const MPM_KRYLOV_VECTOR<TV>& dv=debug_cast<const MPM_KRYLOV_VECTOR<TV>&>(Bdv);
    tmp1=dv;
    if(h) system.forced_collisions.Clean_Memory();
    Make_Feasible(tmp1);
    Compute_Unconstrained(tmp1,h,g,e);
    if(g) Project_Gradient_And_Prune_Constraints(*g,h);
    if(h){
        for(int i=0;i<system.collisions.m;i++){
            const COLLISION& c=system.collisions(i);
            system.forced_collisions.Insert(c.p,c.object);}}

//    if(g) system.Sanity(*g,"g");
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Compute_Unconstrained(const KRYLOV_VECTOR_BASE<T>& Bdv,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const
{
    bool use_midpoint=system.example.use_midpoint;
    T dt=system.example.dt;
    
    const MPM_KRYLOV_VECTOR<TV>& dv=debug_cast<const MPM_KRYLOV_VECTOR<TV>&>(Bdv);
    v1.Copy(use_midpoint?(T).5:1,dv,v0);
    Update_F(v1);
    T midpoint_factor=use_midpoint?(T).25:1;

    system.example.Precompute_Forces(system.example.time,system.example.dt,h);
    if(e){
        T energy=midpoint_factor*system.Inner_Product(dv,dv)/2;
        energy+=system.example.Potential_Energy(system.example.time);
        *e=energy;}

    if(g){
        MPM_KRYLOV_VECTOR<TV>& gg=debug_cast<MPM_KRYLOV_VECTOR<TV>&>(*g);
        tmp2*=0;
        system.example.Add_Forces(tmp2.u,system.example.time);
#pragma omp parallel for
        for(int i=0;i<system.example.valid_grid_indices.m;i++){
            int p=system.example.valid_grid_indices(i);
            gg.u.array(p)=(dv.u.array(p)-dt/system.example.mass.array(p)*tmp2.u.array(p))*midpoint_factor;}}
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Update_F(const MPM_KRYLOV_VECTOR<TV>& v) const
{
    struct HELPER
    {
        MATRIX<T,TV::m> grad_Vp;
        TV Vp,Vn_interpolate;
    };

    //if(system.example.particles.store_S)
    system.example.force_helper.B.Resize(system.example.particles.number);
    system.example.gather_scatter.template Gather<HELPER>(true,
        [](int p,HELPER& h)
        {
            h.grad_Vp=MATRIX<T,TV::m>();
            h.Vp=TV();
            h.Vn_interpolate=TV();
        },
        [this,&v](int p,const PARTICLE_GRID_ITERATOR<TV>& it,HELPER& h)
        {
            int i=v.u.Standard_Index(it.Index());
            h.Vp+=it.Weight()*v.u.array(i);
            h.Vn_interpolate+=it.Weight()*v0.u.array(i);
            h.grad_Vp+=MATRIX<T,TV::m>::Outer_Product(v.u.array(i),it.Gradient());
        },
        [this](int p,HELPER& h)
        {
            MATRIX<T,TV::m> B=system.example.dt*h.grad_Vp,A=B+1;
            system.example.force_helper.B(p)=B; //new addition
            if(system.example.quad_F_coeff) A+=sqr(B)*system.example.quad_F_coeff;
            if(system.example.use_midpoint) A=(A+1)/2;
            system.example.particles.F(p)=A*F0(p);
            if(system.example.particles.store_S){
                //system.example.force_helper.B(p)=B;
                system.example.particles.S(p)=(SYMMETRIC_MATRIX<T,TV::m>::Conjugate(A,S0(p))+system.example.dt*system.example.inv_Wi)/(1+system.example.dt*system.example.inv_Wi);}
            if(system.example.use_midpoint) system.example.particles.X(p)=X0(p)+system.example.dt/2*(h.Vp+h.Vn_interpolate);
            else system.example.particles.X(p)=X0(p)+system.example.dt*h.Vp;
        });
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Restore_F() const
{
#pragma omp parallel for
    for(int k=0;k<system.example.simulated_particles.m;k++){
        int p=system.example.simulated_particles(k);
        system.example.particles.F(p)=F0(p);
        if(system.example.particles.store_S) system.example.particles.S(p)=S0(p);
        system.example.particles.X(p)=X0(p);}
}
//#####################################################################
// Function Adjust_For_Collision
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Adjust_For_Collision(KRYLOV_VECTOR_BASE<T>& Bdv) const
{
    if(!system.example.collision_objects.m) return;
    MPM_KRYLOV_VECTOR<TV>& dv=debug_cast<MPM_KRYLOV_VECTOR<TV>&>(Bdv);
    system.collisions.Remove_All();
    T midpoint_scale=system.example.use_midpoint?(T).5:1;
    system.stuck_nodes.Remove_All();
    system.stuck_velocity.Remove_All();
    T t0=system.example.time;
    T t1=t0+system.example.dt;

#pragma omp parallel for
    for(int tid=0;tid<system.example.gather_scatter.threads;tid++){
        int a=tid*system.example.valid_grid_indices.m/system.example.gather_scatter.threads;
        int b=(tid+1)*system.example.valid_grid_indices.m/system.example.gather_scatter.threads;
        ARRAY<COLLISION> thread_collisions;
        ARRAY<int> thread_stuck_nodes;
        ARRAY<TV> thread_stuck_velocity;
        for(int k=a;k<b;k++){
            int i=system.example.valid_grid_indices(k);
            TV V=v0.u.array(i)+midpoint_scale*dv.u.array(i);
            TV X0=system.example.location.array(i);
            TV X=X0+system.example.dt*V;
            T deepest_phi=FLT_MAX;
            int deepest_index=-1;
            MPM_COLLISION_OBJECT<TV>* deepest_io=0;
            bool stuck=false;
            if(int* object_index=system.forced_collisions.Get_Pointer(i)){
                deepest_index=*object_index;
                deepest_io=system.example.collision_objects(deepest_index);
                T phi0=deepest_io->Phi(X0,t0),phi=deepest_io->Phi(X,t1)-min(phi0,(T)0);
                deepest_phi=phi;
                COLLISION_TYPE type=system.example.collision_objects(deepest_index)->type;
                if(type==COLLISION_TYPE::stick) stuck=true;

            bool allow_sep=system.example.collision_objects(deepest_index)->type==COLLISION_TYPE::separate;
            COLLISION c={deepest_index,i,deepest_phi,0,deepest_io->Normal(X,t1),TV(),deepest_io->Hessian(X,t1),allow_sep};
            thread_collisions.Append(c);
            X-=deepest_phi*c.n;

}
            for(int j=0;j<system.example.collision_objects.m && !stuck;j++){
                MPM_COLLISION_OBJECT<TV>* io=system.example.collision_objects(j);
                T phi0=io->Phi(X0,t0),phi=io->Phi(X,t1)-min(phi0,(T)0);
                COLLISION_TYPE type=system.example.collision_objects(j)->type;
                if(type==COLLISION_TYPE::stick){if(phi0>0) continue;stuck=true;}
                if(type==COLLISION_TYPE::slip && phi0>0) continue;
                if(type==COLLISION_TYPE::separate && phi>collision_thickness) continue;
                if(type!=COLLISION_TYPE::stick && phi>FLT_MAX) continue;
                deepest_index=j;
                deepest_io=io;
                deepest_phi=phi;

            bool allow_sep=system.example.collision_objects(deepest_index)->type==COLLISION_TYPE::separate;
            COLLISION c={deepest_index,i,deepest_phi,0,deepest_io->Normal(X,t1),TV(),deepest_io->Hessian(X,t1),allow_sep};
            thread_collisions.Append(c);
            X-=deepest_phi*c.n;
}
            if(deepest_index==-1) continue;
            if(stuck){
                TV SV=deepest_io->Velocity(X,t1);
                thread_stuck_nodes.Append(i);
                thread_stuck_velocity.Append(SV);
                dv.u.array(i)=(SV-v0.u.array(i))/midpoint_scale; // Come to rest
                continue;}




            V=(X-system.example.location.array(i))/system.example.dt;
            dv.u.array(i)=(V-v0.u.array(i))/midpoint_scale;}
#pragma omp critical
        {
            if(!system.stuck_nodes.m) system.stuck_nodes.Exchange(thread_stuck_nodes);
            else system.stuck_nodes.Append_Elements(thread_stuck_nodes);
            if(!system.stuck_velocity.m) system.stuck_velocity.Exchange(thread_stuck_velocity);
            else system.stuck_velocity.Append_Elements(thread_stuck_velocity);
            if(!system.collisions.m) system.collisions.Exchange(thread_collisions);
            else system.collisions.Append_Elements(thread_collisions);
        }
    }
}
//#####################################################################
// Function Project_Gradient_And_Prune_Constraints
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Project_Gradient_And_Prune_Constraints(KRYLOV_VECTOR_BASE<T>& Bg,bool allow_sep) const
{
    if(!system.collisions.m && !system.stuck_nodes.m) return;
    MPM_KRYLOV_VECTOR<TV>& g=debug_cast<MPM_KRYLOV_VECTOR<TV>&>(Bg);
#pragma omp parallel for
    for(int i=0;i<system.stuck_nodes.m;i++)
        g.u.array(system.stuck_nodes(i))=TV();

#pragma omp parallel for
    for(int i=0;i<system.collisions.m;i++){
        COLLISION& c=system.collisions(i);
        TV &gv=g.u.array(c.p);
        c.n_dE=gv.Dot(c.n);
        if(!allow_sep || !c.allow_sep || c.n_dE>=0){
            c.H_dE=c.H*gv;
            gv-=c.n_dE*c.n+c.phi*c.H_dE;}}

    for(int i=system.collisions.m-1;i>=0;i--){
        const COLLISION& c=system.collisions(i);
        if(allow_sep && c.allow_sep && c.n_dE<0)
            system.collisions.Remove_Index_Lazy(i);}

    LOG::printf("collisions: %i\n",system.collisions.m);
}
//#####################################################################
// Function Make_Feasible
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Make_Feasible(KRYLOV_VECTOR_BASE<T>& dv) const
{
//    system.Sanity(dv,"before make feasible");

    tmp2.Copy(0,tmp2,dv);
    Adjust_For_Collision(dv);
    tmp2.Copy(-1,dv,tmp2);
    system.Sanity(tmp2,"make feasible diff");

//    system.Sanity(dv,"after make feasible");
}
//#####################################################################
// Function Initial_Guess
//#####################################################################
template<class TV> bool MPM_OBJECTIVE<TV>::
Initial_Guess(KRYLOV_VECTOR_BASE<T>& Bdv,T tolerance,bool no_test) const
{
    if(no_test) LOG::printf("NO TEST INITIAL GUESS\n");
    MPM_KRYLOV_VECTOR<TV>& dv=debug_cast<MPM_KRYLOV_VECTOR<TV>&>(Bdv);
    system.forced_collisions.Remove_All();
    T e0=0,e1=0;
    dv*=0;
    Compute(dv,0,&tmp0,&e0);
    T factor=system.example.use_midpoint?4:1;

    dv.Copy(-factor,tmp0);
    if(no_test) return true;
    T norm_grad0=sqrt(system.Inner_Product(tmp0,tmp0));
    Compute(dv,0,&tmp0,&e1);
    T norm_grad1=sqrt(system.Inner_Product(tmp0,tmp0));
    if(e1<max(e0,e0*(T)1.1) || norm_grad1<norm_grad0*(T)1.1 || norm_grad1<tolerance) return true;
    dv*=0;
    return false;
}
//#####################################################################
// Function Reset
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Reset()
{
    F0.Resize(system.example.particles.X.m);
    if(system.example.particles.store_S) S0.Resize(system.example.particles.X.m);
    X0.Resize(system.example.particles.X.m);
#pragma omp parallel for
    for(int k=0;k<system.example.simulated_particles.m;k++){
        int p=system.example.simulated_particles(k);
        F0(p)=system.example.particles.F(p);
        if(system.example.particles.store_S) S0(p)=system.example.particles.S(p);
        X0(p)=system.example.particles.X(p);}

    v0.u=system.example.velocity;
    v1.u.Resize(v0.u.domain);
    tmp0.u.Resize(v0.u.domain);
    tmp1.u.Resize(v0.u.domain);
    tmp2.u.Resize(v0.u.domain);
}
//#####################################################################
// Function Test_Diff
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Test_Diff(const KRYLOV_VECTOR_BASE<T>& dv)
{
    KRYLOV_VECTOR_BASE<T>& ddv=*dv.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& g0=*dv.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& g1=*dv.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& a=*dv.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& b=*dv.Clone_Default();
    KRYLOV_VECTOR_BASE<T>& t0=*dv.Clone_Default();

    t0=dv;
    this->Test(t0,system);

    T eps=(T)1e-6;
    static RANDOM_NUMBERS<T> random;

    T e0=0,e1=0;
    Compute(t0,&system,&g0,&e0);

    for(int i=0,n=ddv.Raw_Size();i<n;i++)
        ddv.Raw_Get(i)=random.Get_Uniform_Number(-eps,eps);
    system.Project(ddv);
    b.Copy(1,t0,ddv);

    ddv.Copy(-1,t0,b);

    system.Multiply(ddv,a);
    system.Project(a);

    // new dv
    Compute(b,&system,&g1,&e1);
    system.Multiply(ddv,b);
    system.Project(b);

    T test0a=(system.Inner_Product(g1,ddv)+system.Inner_Product(g0,ddv))/(2*eps);
    T test0b=(e1-e0)/eps;
    T test0=(test0b-test0a)/max(abs(test0a),(T)1e-20);

    a.Copy(1,b,a);
    a.Copy((T).5,a);
    T test1a=sqrt(system.Inner_Product(a,a))/eps;

    b.Copy(-1,g0,g1);
    T test1b=sqrt(system.Inner_Product(b,b))/eps;

    a.Copy(-1,b,a);
    T test1=sqrt(system.Inner_Product(a,a))/(eps*max(abs(test1a),(T)1e-20));

    if((e0!=0)&&(e1-e0!=0)) LOG::cout<<"RELATIVE ENERGY "<<-log10(abs(e1-e0)/abs(e0))<<std::endl;
    LOG::cout<<"energy diff test C "<<test0a<<"    "<<test0b<<"    "<<test0<<std::endl;
    LOG::cout<<"force diff test C "<<test1a<<"    "<<test1b<<"    "<<test1<<std::endl;

    delete &ddv;
    delete &g0;
    delete &g1;
    delete &a;
    delete &b;
    delete &t0;
}
template class MPM_OBJECTIVE<VECTOR<float,2> >;
template class MPM_OBJECTIVE<VECTOR<float,3> >;
template class MPM_OBJECTIVE<VECTOR<double,2> >;
template class MPM_OBJECTIVE<VECTOR<double,3> >;
}
