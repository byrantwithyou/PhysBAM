//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
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

    system.example.Precompute_Forces(system.example.time);
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

    system.example.gather_scatter.template Gather<HELPER>(
        [](int p,HELPER& h)
        {
            h.grad_Vp=MATRIX<T,TV::m>();
            h.Vp=TV();
            h.Vn_interpolate=TV();
        },
        [this,&v](int p,const PARTICLE_GRID_ITERATOR<TV>& it,HELPER& h)
        {
            h.Vp+=it.Weight()*v.u(it.Index());
            h.Vn_interpolate+=it.Weight()*v0.u(it.Index());
            h.grad_Vp+=MATRIX<T,TV::m>::Outer_Product(v.u(it.Index()),it.Gradient());
        },
        [this](int p,HELPER& h)
        {
            system.example.particles.F(p)=F0(p)+system.example.dt/(system.example.use_midpoint+1)*h.grad_Vp*F0(p);
            if(system.example.use_midpoint) system.example.particles.X(p)=X0(p)+system.example.dt/2*(h.Vp+h.Vn_interpolate);
            else system.example.particles.X(p)=X0(p)+system.example.dt*h.Vp;
        },true);
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

    // TODO: parallelize
    for(int k=0;k<system.example.valid_grid_indices.m;k++){
        int i=system.example.valid_grid_indices(k);
        TV V=v0.u.array(i)+midpoint_scale*dv.u.array(i);
        TV X0=system.example.location.array(i);
        TV X=X0+system.example.dt*V;
        T deepest_phi=collision_thickness;
        int deepest_index=-1;
        IMPLICIT_OBJECT<TV>* deepest_io=0;
        if(int* object_index=system.forced_collisions.Get_Pointer(i)){
            deepest_index=*object_index;
            deepest_io=system.example.collision_objects(deepest_index).io;
            T phi0=deepest_io->Extended_Phi(X0),phi=deepest_io->Extended_Phi(X);
            if(phi0<0) phi-=phi0;
            deepest_phi=phi;}
        for(int j=0;j<system.example.collision_objects.m;j++){
            IMPLICIT_OBJECT<TV>* io=system.example.collision_objects(j).io;
            T phi0=io->Extended_Phi(X0),phi=io->Extended_Phi(X);
            if(phi0<0) phi-=phi0;
            if(phi<deepest_phi){
                deepest_phi=phi;
                deepest_io=io;
                deepest_index=j;}}
        if(deepest_index==-1) continue;
        if(system.example.collision_objects(deepest_index).sticky){
            system.stuck_nodes.Append(i);
            continue;}
        COLLISION c={deepest_index,i,deepest_phi,0,deepest_io->Extended_Normal(X),TV(),deepest_io->Hessian(X)};
        system.collisions.Append(c);
        X-=deepest_phi*c.n;
        V=(X-system.example.location.array(i))/system.example.dt;
        dv.u.array(i)=(V-v0.u.array(i))/midpoint_scale;}
}
//#####################################################################
// Function Project_Gradient_And_Prune_Constraints
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Project_Gradient_And_Prune_Constraints(KRYLOV_VECTOR_BASE<T>& Bg,bool allow_sep) const
{
    if(!system.collisions.m && !system.stuck_nodes.m) return;
    MPM_KRYLOV_VECTOR<TV>& g=debug_cast<MPM_KRYLOV_VECTOR<TV>&>(Bg);
    g.u.array.Subset(system.stuck_nodes).Fill(TV());

#pragma omp parallel for
    for(int i=0;i<system.collisions.m;i++){
        COLLISION& c=system.collisions(i);
        TV &gv=g.u.array(c.p);
        c.n_dE=gv.Dot(c.n);
        if(!allow_sep || c.n_dE>=0){
            c.H_dE=c.H*gv;
            gv-=c.n_dE*c.n+c.phi*c.H_dE;}}

    for(int i=system.collisions.m-1;i>=0;i--){
        COLLISION& c=system.collisions(i);
        if(allow_sep && c.n_dE<0)
            system.collisions.Remove_Index_Lazy(i);}
}
//#####################################################################
// Function Make_Feasible
//#####################################################################
template<class TV> void MPM_OBJECTIVE<TV>::
Make_Feasible(KRYLOV_VECTOR_BASE<T>& dv) const
{
    Adjust_For_Collision(dv);
}
//#####################################################################
// Function Initial_Guess
//#####################################################################
template<class TV> bool MPM_OBJECTIVE<TV>::
Initial_Guess(KRYLOV_VECTOR_BASE<T>& Bdv,T tolerance) const
{
    MPM_KRYLOV_VECTOR<TV>& dv=debug_cast<MPM_KRYLOV_VECTOR<TV>&>(Bdv);
    T e0=0,e1=0;
    dv*=0;
    Compute(dv,0,&tmp0,&e0);
    T factor=system.example.use_midpoint?4:1;

    dv.Copy(-factor,tmp0);
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
    X0.Resize(system.example.particles.X.m);
#pragma omp parallel for
    for(int k=0;k<system.example.simulated_particles.m;k++){
        int p=system.example.simulated_particles(k);
        F0(p)=system.example.particles.F(p);
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
