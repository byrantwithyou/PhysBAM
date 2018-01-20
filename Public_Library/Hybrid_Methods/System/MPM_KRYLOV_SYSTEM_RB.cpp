//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE_RB.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_SYSTEM_RB.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR_RB.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_KRYLOV_SYSTEM_RB<TV>::
MPM_KRYLOV_SYSTEM_RB(MPM_EXAMPLE_RB<TV>& example)
    :KRYLOV_SYSTEM_BASE<T>(false,false),example(example),tmp(*new MPM_KRYLOV_VECTOR_RB<TV>(example.valid_grid_indices))
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_KRYLOV_SYSTEM_RB<TV>::
~MPM_KRYLOV_SYSTEM_RB()
{
    delete &tmp;
}
//#####################################################################
// Function Sanity
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM_RB<TV>::
Sanity(const KRYLOV_VECTOR_BASE<T>& v,const char* str) const
{
    return;
    const MPM_KRYLOV_VECTOR_RB<TV>& vv=debug_cast<const MPM_KRYLOV_VECTOR_RB<TV>&>(v);
    LOG::printf("%s: ",str);
    for(int i=0;i<collisions.m;i++){
        auto& c=collisions(i);
        LOG::printf("%P ",vv.u.array(c.p).Dot(c.n));}
    LOG::printf("\n");
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM_RB<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF,bool transpose) const
{
    const MPM_KRYLOV_VECTOR_RB<TV>& V=debug_cast<const MPM_KRYLOV_VECTOR_RB<TV>&>(BV);
    MPM_KRYLOV_VECTOR_RB<TV>& F=debug_cast<MPM_KRYLOV_VECTOR_RB<TV>&>(BF);
    tmp=V;
#pragma omp parallel for
    for(int i=0;i<stuck_nodes.m;i++)
        tmp.u.array(stuck_nodes(i))=TV();
#pragma omp parallel for
    for(int i=0;i<collisions.m;i++){
        const COLLISION& c=collisions(i);
        tmp.u.array(c.p).Project_Orthogonal_To_Unit_Direction(c.n);}
    
//    LOG::printf("project: %i\n",collisions.m);

    F*=0;
    example.Add_Hessian_Times(F.u,tmp.u,F.twists,tmp.twists,example.time,transpose);
    example.Reflection_Boundary_Condition(F.u,true);
    T scale=example.use_midpoint?(T).25:1,scaled_dt_squared=sqr(example.dt*scale);

#pragma omp parallel for
    for(int k=0;k<example.valid_grid_indices.m;k++){
        int i=example.valid_grid_indices(k);
        F.u.array(i)=scaled_dt_squared/example.mass.array(i)*F.u.array(i)+tmp.u.array(i)*scale;}

#pragma omp parallel for
    for(int i=0;i<F.twists.m;i++){
        F.twists(i)=rigid_mass_inverse(i)*F.twists(i)*scaled_dt_squared+tmp.twists(i)*scale;}

#pragma omp parallel for
    for(int i=0;i<stuck_nodes.m;i++)
        F.u.array(stuck_nodes(i))=TV();
#pragma omp parallel for
    for(int i=0;i<collisions.m;i++){
        const COLLISION& c=collisions(i);
        TV& v=F.u.array(c.p),tt=V.u.array(c.p); // tmp -> V?
        v-=c.n*c.n.Dot(v)+example.dt*(c.H*tt*c.n_dE+c.n.Dot(tt)*c.H_dE+c.H_dE.Dot(tt)*c.n);}

//    Sanity(BF,"H");
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MPM_KRYLOV_SYSTEM_RB<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const MPM_KRYLOV_VECTOR_RB<TV>& X=debug_cast<const MPM_KRYLOV_VECTOR_RB<TV>&>(x);
    const MPM_KRYLOV_VECTOR_RB<TV>& Y=debug_cast<const MPM_KRYLOV_VECTOR_RB<TV>&>(y);
    T r=0;
#pragma omp parallel for reduction(+:r)
    for(int t=0;t<example.threads;t++){
        int a=t*example.valid_grid_indices.m/example.threads;
        int b=(t+1)*example.valid_grid_indices.m/example.threads;
        for(int k=a;k<b;k++){
            int i=example.valid_grid_indices(k);
            r+=example.mass.array(i)*X.u.array(i).Dot(Y.u.array(i));}}
#pragma omp parallel for reduction(+:r)
    for(int t=0;t<example.threads;t++){
        int a=t*X.twists.m/example.threads;
        int b=(t+1)*X.twists.m/example.threads;
        for(int k=a;k<b;k++){
            r+=rigid_mass(k).Inner_Product(X.twists(k),Y.twists(k));}}
    return r;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MPM_KRYLOV_SYSTEM_RB<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    return sqrt(Inner_Product(BR,BR));
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM_RB<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM_RB<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM_RB<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BR) const
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM_RB<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
}
template class MPM_KRYLOV_SYSTEM_RB<VECTOR<float,2> >;
template class MPM_KRYLOV_SYSTEM_RB<VECTOR<float,3> >;
template class MPM_KRYLOV_SYSTEM_RB<VECTOR<double,2> >;
template class MPM_KRYLOV_SYSTEM_RB<VECTOR<double,3> >;
}
