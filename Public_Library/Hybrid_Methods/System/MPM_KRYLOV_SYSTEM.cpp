//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_SYSTEM.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_VECTOR.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_KRYLOV_SYSTEM<TV>::
MPM_KRYLOV_SYSTEM(MPM_EXAMPLE<TV>& example)
    :KRYLOV_SYSTEM_BASE<T>(false,false),example(example),tmp(*new MPM_KRYLOV_VECTOR<TV>(example.valid_grid_indices))
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_KRYLOV_SYSTEM<TV>::
~MPM_KRYLOV_SYSTEM()
{
    delete &tmp;
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const MPM_KRYLOV_VECTOR<TV>& V=debug_cast<const MPM_KRYLOV_VECTOR<TV>&>(BV);
    MPM_KRYLOV_VECTOR<TV>& F=debug_cast<MPM_KRYLOV_VECTOR<TV>&>(BF);
    tmp=V;
    tmp.u.array.Subset(stuck_nodes).Fill(TV());
    for(int i=0;i<collisions.m;i++){
        const COLLISION& c=collisions(i);
        tmp.u.array(c.p).Project_Orthogonal_To_Unit_Direction(c.n);}

    F*=0;
    example.Add_Hessian_Times(F.u,tmp.u,example.time);
    T scale=example.use_midpoint?(T).25:1,scaled_dt_squared=sqr(example.dt*scale);

#pragma omp parallel for
    for(int k=0;k<example.valid_grid_indices.m;k++){
        int i=example.valid_grid_indices(k);
        F.u.array(i)=scaled_dt_squared/example.mass.array(i)*F.u.array(i)+tmp.u.array(i)*scale;}

    F.u.array.Subset(stuck_nodes).Fill(TV());
    for(int i=0;i<collisions.m;i++){
        const COLLISION& c=collisions(i);
        TV& v=F.u.array(c.p),tt=V.u.array(c.p); // tmp -> V?
        v-=c.n*c.n.Dot(v)+example.dt*(c.H*tt*c.n_dE+c.n.Dot(tt)*c.H_dE+c.H_dE.Dot(tt)*c.n);}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MPM_KRYLOV_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const MPM_KRYLOV_VECTOR<TV>& X=debug_cast<const MPM_KRYLOV_VECTOR<TV>&>(x);
    const MPM_KRYLOV_VECTOR<TV>& Y=debug_cast<const MPM_KRYLOV_VECTOR<TV>&>(y);
    T r=0;
#pragma omp parallel for reduction(+:r)
    for(int t=0;t<example.threads;t++){
        int a=t*example.valid_grid_indices.m/example.threads;
        int b=(t+1)*example.valid_grid_indices.m/example.threads;
        for(int k=a;k<b;k++){
            int i=example.valid_grid_indices(k);
            r+=example.mass.array(i)*X.u.array(i).Dot(Y.u.array(i));}}
    return r;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MPM_KRYLOV_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    return sqrt(Inner_Product(BR,BR));
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BR) const
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void MPM_KRYLOV_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
}
template class MPM_KRYLOV_SYSTEM<VECTOR<float,2> >;
template class MPM_KRYLOV_SYSTEM<VECTOR<float,3> >;
template class MPM_KRYLOV_SYSTEM<VECTOR<double,2> >;
template class MPM_KRYLOV_SYSTEM<VECTOR<double,3> >;
}
