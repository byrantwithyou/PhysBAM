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
#include <Hybrid_Methods/Forces/MPM_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_FINITE_ELEMENTS<TV>::
MPM_FINITE_ELEMENTS(const MPM_FORCE_HELPER<TV>& force_helper,
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model,
    GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles)
    :BASE(force_helper),constitutive_model(constitutive_model),affect_all(!affected_particles),
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
template<class TV> MPM_FINITE_ELEMENTS<TV>:: 
~MPM_FINITE_ELEMENTS()
{
    if(!affect_all){
        delete &gather_scatter.simulated_particles;
        delete &gather_scatter;}
    delete &constitutive_model;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void MPM_FINITE_ELEMENTS<TV>:: 
Precompute(const T time,const T dt,bool want_dE,bool want_ddE)
{
    PHYSBAM_ASSERT(!want_ddE||want_dE);
    U.Resize(particles.X.m);
    FV.Resize(particles.X.m);
    sigma.Resize(particles.X.m);
    dPi_dF.Resize(particles.X.m);
    PFT.Resize(particles.X.m);

#pragma omp parallel for
    for(int k=0;k<gather_scatter.simulated_particles.m;k++){
        int p=gather_scatter.simulated_particles(k);
        MATRIX<T,TV::m> V_local;
        particles.F(p).Fast_Singular_Value_Decomposition(U(p),sigma(p),V_local);
        FV(p)=force_helper.Fn(p)*V_local;
        PFT(p)=U(p)*constitutive_model.P_From_Strain(sigma(p),p).Times_Transpose(FV(p));
        constitutive_model.Isotropic_Stress_Derivative(sigma(p),dPi_dF(p),p);}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_FINITE_ELEMENTS<TV>:: 
Potential_Energy(const T time) const
{
    T pe=0;
#pragma omp parallel for reduction(+:pe)
    for(int k=0;k<gather_scatter.simulated_particles.m;k++){
        int p=gather_scatter.simulated_particles(k);
        pe+=particles.volume(p)*constitutive_model.Energy_Density(sigma(p),p);}
    return pe;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_FINITE_ELEMENTS<TV>:: 
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    T c=force_helper.quad_F_coeff;
    bool use_c=c!=0;
    gather_scatter.template Scatter<MATRIX<T,TV::m> >(
        [this,c,use_c](int p,MATRIX<T,TV::m>& AF)
        {
            MATRIX<T,TV::m> P=PFT(p);
            if(use_c){MATRIX<T,TV::m> B=force_helper.B(p);P+=(P.Times_Transpose(B)+B.Transpose_Times(P))*c;}
            AF=P*particles.volume(p);
        },
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,const MATRIX<T,TV::m>& A)
        {F(it.Index())-=A*it.Gradient();},
        [](int p,const MATRIX<T,TV::m>& A){},true);
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_FINITE_ELEMENTS<TV>:: 
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& Z,const T time) const
{
    tmp.Resize(particles.X.m);
    
    T c=force_helper.quad_F_coeff;
    bool use_c=c!=0;
    gather_scatter.template Gather<int>(
        [this](int p,int data){tmp(p)=MATRIX<T,TV::m>();},
        [this,&Z](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {
            tmp(p)+=MATRIX<T,TV::m>::Outer_Product(Z(it.Index()),it.Gradient());
        },
        [this,c,use_c](int p,int data){
            MATRIX<T,TV::m> W=tmp(p);
            if(use_c){
                MATRIX<T,TV::m> B=force_helper.B(p);
                W+=(W*B+B*W)*c;}
            MATRIX<T,TV::m> UU=U(p),FVV=FV(p),G=UU.Transpose_Times(W*FVV);
            MATRIX<T,TV::m> M=UU*dPi_dF(p).Differential(G).Times_Transpose(FVV);
            if(use_c){
                MATRIX<T,TV::m> B=force_helper.B(p),P=PFT(p);
                M+=(M.Times_Transpose(B)+B.Transpose_Times(M)+P.Times_Transpose(tmp(p))+tmp(p).Transpose_Times(P))*c;}
            tmp(p)=M*particles.volume(p);
        },true);

    gather_scatter.template Scatter<int>(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {F(it.Index())+=tmp(p)*it.Gradient();},true);
}
template class MPM_FINITE_ELEMENTS<VECTOR<float,1> >;
template class MPM_FINITE_ELEMENTS<VECTOR<float,2> >;
template class MPM_FINITE_ELEMENTS<VECTOR<float,3> >;
template class MPM_FINITE_ELEMENTS<VECTOR<double,1> >;
template class MPM_FINITE_ELEMENTS<VECTOR<double,2> >;
template class MPM_FINITE_ELEMENTS<VECTOR<double,3> >;
}
