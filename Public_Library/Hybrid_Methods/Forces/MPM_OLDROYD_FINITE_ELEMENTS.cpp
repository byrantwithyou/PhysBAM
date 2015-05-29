//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_OLDROYD_FINITE_ELEMENTS<TV>::
MPM_OLDROYD_FINITE_ELEMENTS(MPM_PARTICLES<TV>& particles,OLDROYD_CONSTITUTIVE_MODEL<TV>& constitutive_model,
    GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles)
    :BASE(particles),constitutive_model(constitutive_model),affect_all(!affected_particles),
    gather_scatter(affect_all?gather_scatter_input:*new GATHER_SCATTER<TV>(gather_scatter_input.grid,*new ARRAY<int>(*affected_particles)))
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
Capture_Stress()
{
    Fn=particles.F;
    Sn=particles.S;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void MPM_OLDROYD_FINITE_ELEMENTS<TV>:: 
Precompute(const T time)
{
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
    gather_scatter.template Scatter<MATRIX<T,TV::m> >(
        [this](int p,MATRIX<T,TV::m>& A)
        {
            SYMMETRIC_MATRIX<T,TV::m> dQ;
            MATRIX<T,TV::m> dP;
            constitutive_model.Gradient(dP,dQ,p);
            A=(dP.Times_Transpose(Fn(p))+dQ.Times_Transpose(Sn(p)*2))*particles.volume(p);
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

    gather_scatter.template Gather<int>(
        [this](int p,int data){tmp(p)=MATRIX<T,TV::m>();},
        [this,&V](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {tmp(p)+=MATRIX<T,TV::m>::Outer_Product(V(it.Index()),it.Gradient());},
        [this](int p,int data){
            SYMMETRIC_MATRIX<T,TV::m> B=(tmp(p)*Sn(p)).Twice_Symmetric_Part(),dQ;
            MATRIX<T,TV::m> C=tmp(p)*Fn(p),dP;
            constitutive_model.Hessian(C,B,dP,dQ,p);
            tmp(p)=(dP.Times_Transpose(Fn(p))+dQ.Times_Transpose(Sn(p)*2))*particles.volume(p);},true);

    gather_scatter.template Scatter<int>(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {F(it.Index())+=tmp(p)*it.Gradient();},true);
}
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<float,2> >;
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<float,3> >;
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<double,2> >;
template class MPM_OLDROYD_FINITE_ELEMENTS<VECTOR<double,3> >;
}
