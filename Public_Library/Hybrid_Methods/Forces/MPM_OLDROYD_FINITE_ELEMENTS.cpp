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
    GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles,const T& inv_Wi,const T viscosity)
    :BASE(force_helper),constitutive_model(constitutive_model),affect_all(!affected_particles),
    gather_scatter(affect_all?gather_scatter_input:*new GATHER_SCATTER<TV>(gather_scatter_input.grid,*new ARRAY<int>(*affected_particles))),stored_dt(0),inv_Wi(inv_Wi),viscosity(viscosity)
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
Precompute(const T time,const T dt)
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
    T q_scale=(T).5/sqr(stored_dt);
#pragma omp parallel for reduction(+:pe)
    for(int k=0;k<gather_scatter.simulated_particles.m;k++){
        int p=gather_scatter.simulated_particles(k);
        pe+=particles.volume(p)*constitutive_model.Energy_Density(p);
        if(viscosity){
            MATRIX<T,TV::m> A=force_helper.A(p)-MATRIX<T,TV::m>::Identity_Matrix();
            T q=(A.Frobenius_Norm_Squared()+MATRIX<T,TV::m>::Inner_Product(A.Transposed(),A))*q_scale;
            T c=force_helper.Fn(p).Determinant()*viscosity*q;
            pe+=particles.volume(p)*c;}}
    return pe;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_OLDROYD_FINITE_ELEMENTS<TV>:: 
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    T zeta=2/(1+stored_dt*inv_Wi);
    T q_scale=(T)1/sqr(stored_dt);
    gather_scatter.template Scatter<MATRIX<T,TV::m> >(
        [this,zeta,q_scale](int p,MATRIX<T,TV::m>& A)
        {
            SYMMETRIC_MATRIX<T,TV::m> dQ;
            MATRIX<T,TV::m> dP;
            constitutive_model.Gradient(dP,dQ,p);
            A=(dP.Times_Transpose(force_helper.Fn(p))+zeta*(dQ*force_helper.A(p)*force_helper.Sn(p)))*particles.volume(p);
            if(viscosity){
                MATRIX<T,TV::m> Ap=force_helper.A(p)-MATRIX<T,TV::m>::Identity_Matrix();
                A+=(Ap+Ap.Transposed())*(q_scale*viscosity*force_helper.Fn(p).Determinant()*particles.volume(p));}
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
    tmp2.Resize(particles.X.m);

    T zeta=2/(1+stored_dt*inv_Wi);
    T q_scale=(T)1/sqr(stored_dt);
    gather_scatter.template Gather<int>(
        [this](int p,int data){tmp(p)=MATRIX<T,TV::m>();tmp2(p)=MATRIX<T,TV::m>();},
        [this,&V](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {
            tmp(p)+=MATRIX<T,TV::m>::Outer_Product(V(it.Index()),it.Gradient());
            if(viscosity){
                tmp2(p)+=MATRIX<T,TV::m>::Outer_Product(V(it.Index()),it.Gradient());
                tmp2(p)+=MATRIX<T,TV::m>::Outer_Product(it.Gradient(),V(it.Index()));}
        },
        [this,zeta,q_scale](int p,int data){
            SYMMETRIC_MATRIX<T,TV::m> B=(tmp(p)*force_helper.Sn(p)).Twice_Symmetric_Part(),dQ2;
            MATRIX<T,TV::m> dQ;
            MATRIX<T,TV::m> C=tmp(p)*force_helper.Fn(p),dP,dP2;
            constitutive_model.Hessian(C,B,dP,dQ,p);
            constitutive_model.Gradient(dP2,dQ2,p);
            MATRIX<T,TV::m> G=dP.Times_Transpose(force_helper.Fn(p));
            //G+=zeta*dQ*force_helper.Sn(p).Times_Transpose(force_helper.A(p));
            G+=zeta*dQ*force_helper.A(p)*force_helper.Sn(p);
            G+=zeta*dQ2*tmp(p)*force_helper.Sn(p);
            tmp(p)=G*particles.volume(p);
            if(viscosity) tmp2(p)*=(q_scale*viscosity*force_helper.Fn(p).Determinant()*particles.volume(p));},true);

    gather_scatter.template Scatter<int>(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {F(it.Index())+=tmp(p)*it.Gradient()+tmp2(p)*it.Gradient();},true);
}
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<float,2> >;
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<float,3> >;
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<double,2> >;
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<double,3> >;
}
