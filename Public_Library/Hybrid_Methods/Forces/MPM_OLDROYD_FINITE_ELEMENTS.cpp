//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_OLDROYD_FINITE_ELEMENTS<TV>::
MPM_OLDROYD_FINITE_ELEMENTS(MPM_FORCE_HELPER<TV>& force_helper,
    OLDROYD_CONSTITUTIVE_MODEL<TV>& constitutive_model,
    GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles,
    const T& inv_Wi,T quad_F_coeff)
    :BASE(force_helper),constitutive_model(constitutive_model),affect_all(!affected_particles),
    gather_scatter(affect_all?gather_scatter_input:*new GATHER_SCATTER<TV>(gather_scatter_input.grid,*new ARRAY<int>(*affected_particles))),stored_dt(0),inv_Wi(inv_Wi)
{
    if(!affect_all){
        gather_scatter.weights=gather_scatter_input.weights;
        gather_scatter.threads=gather_scatter_input.threads;
        gather_scatter.Prepare_Scatter(particles);}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_OLDROYD_FINITE_ELEMENTS<TV>:: 
~MPM_OLDROYD_FINITE_ELEMENTS()
{
    if(!affect_all){
        delete &gather_scatter.simulated_particles;
        delete &gather_scatter;}
    delete &constitutive_model;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void MPM_OLDROYD_FINITE_ELEMENTS<TV>:: 
Precompute(const T time,const T dt,bool want_dE,bool want_ddE)
{
    PHYSBAM_ASSERT(particles.store_S);
    stored_dt=dt;
    constitutive_model.Resize(particles.number);
#pragma omp parallel for
    for(int k=0;k<gather_scatter.simulated_particles.m;k++){
        int p=gather_scatter.simulated_particles(k);
        constitutive_model.Precompute(particles.F(p),particles.S(p),p);}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_OLDROYD_FINITE_ELEMENTS<TV>:: 
Potential_Energy(const T time) const
{
    T pe=0;
#pragma omp parallel for reduction(+:pe)
    for(int k=0;k<gather_scatter.simulated_particles.m;k++){
        int p=gather_scatter.simulated_particles(k);
        pe+=particles.volume(p)*constitutive_model.Energy_Density(p);}
    return pe;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_OLDROYD_FINITE_ELEMENTS<TV>:: 
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    T zeta=1/(1+stored_dt*inv_Wi);
    T c=force_helper.quad_F_coeff;
    gather_scatter.template Scatter<MATRIX<T,TV::m> >(
        [this,zeta,c](int p,MATRIX<T,TV::m>& AF)
        {
            SYMMETRIC_MATRIX<T,TV::m> dQ,S=force_helper.Sn(p);
            MATRIX<T,TV::m> dP,F=force_helper.Fn(p),B=force_helper.B(p);
            constitutive_model.Gradient(dP,dQ,p);
            MATRIX<T,TV::m> A=B+1+c*sqr(B);
            MATRIX<T,TV::m> N=dP.Times_Transpose(F);
            MATRIX<T,TV::m> Q=2*zeta*dQ*A*S;
            MATRIX<T,TV::m> P=N+Q;
            MATRIX<T,TV::m> TP=(P+c*(P.Times_Transpose(B)+B.Transpose_Times(P)));
            AF=TP*particles.volume(p);
        },
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,const MATRIX<T,TV::m>& A)
        {F(it.Index())-=A*it.Gradient();},
        [](int p,const MATRIX<T,TV::m>& A){},true);
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_OLDROYD_FINITE_ELEMENTS<TV>:: 
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& V,const T time) const
{
    tmp.Resize(particles.X.m);

    T zeta=1/(1+stored_dt*inv_Wi);
    T c=force_helper.quad_F_coeff;
    gather_scatter.template Gather<int>(
        [this](int p,int data){tmp(p)=MATRIX<T,TV::m>();},
        [this,&V](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {
            tmp(p)+=MATRIX<T,TV::m>::Outer_Product(V(it.Index()),it.Gradient());
        },
        [this,zeta,c](int p,int data){
            SYMMETRIC_MATRIX<T,TV::m> dQ,S=force_helper.Sn(p),G;
            MATRIX<T,TV::m> dP,F=force_helper.Fn(p),B=force_helper.B(p),L;
            constitutive_model.Gradient(dP,dQ,p);
            MATRIX<T,TV::m> A=B+1+c*sqr(B);
            MATRIX<T,TV::m> N=dP.Times_Transpose(F);
            MATRIX<T,TV::m> Q=2*zeta*dQ*A*S;
            MATRIX<T,TV::m> P=N+Q;
            MATRIX<T,TV::m> V=(P.Times_Transpose(tmp(p))+tmp(p).Transpose_Times(P))*c;
            MATRIX<T,TV::m> W=tmp(p)+(tmp(p)*B+B*tmp(p))*c;
            MATRIX<T,TV::m> D=W*F;
            SYMMETRIC_MATRIX<T,TV::m> E=zeta*(A*S.Times_Transpose(W)).Twice_Symmetric_Part();
            constitutive_model.Hessian(D,E,L,G,p);
            MATRIX<T,TV::m> H=L.Times_Transpose(F);
            MATRIX<T,TV::m> K=2*zeta*(G*A+dQ*W)*S;
            MATRIX<T,TV::m> M=H+K;
            MATRIX<T,TV::m> Z=M+(M.Times_Transpose(B)+B.Transpose_Times(M))*c;
            tmp(p)=(Z+V)*particles.volume(p);
        },true);

    gather_scatter.template Scatter<int>(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {F(it.Index())+=tmp(p)*it.Gradient();},true);
}
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<float,2> >;
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<float,3> >;
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<double,2> >;
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<double,3> >;
}
