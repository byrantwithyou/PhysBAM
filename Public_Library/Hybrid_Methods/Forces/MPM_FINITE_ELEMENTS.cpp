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
#include <Deformables/Constitutive_Models/MPM_PLASTICITY_MODEL.h>
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
    GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles,
    MPM_PLASTICITY_MODEL<TV>* plasticity)
    :BASE(force_helper),constitutive_model(constitutive_model),plasticity(plasticity),affect_all(!affected_particles),
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
    delete plasticity;
}
//#####################################################################
// Function Transform_Isotropic_Stress_Derivative
//#####################################################################
template<class T> void
Transform_Isotropic_Stress_Derivative(DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,1>& dPi_dF,
        const DIAGONAL_MATRIX<T,1>& P,const DIAGONAL_MATRIX<T,1>& P_hat,
        const MATRIX<T,1>& dstrain,const SYMMETRIC_TENSOR<T,0,1>& ddstrain,
        const MATRIX<T,1,0>& rdstrain,const MATRIX<T,0>& rxstrain,const VECTOR<T,1>& sigma){
    dPi_dF.Set_Hessian_Block(Contract(ddstrain,P.x)+SYMMETRIC_MATRIX<T,1>::Conjugate_With_Transpose(dstrain,dPi_dF.Get_Hessian_Block()));
}
template<class T> void
Transform_Isotropic_Stress_Derivative(DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dPi_dF,
        const DIAGONAL_MATRIX<T,2>& P,const DIAGONAL_MATRIX<T,2>& P_hat,
        const MATRIX<T,2>& dstrain,const SYMMETRIC_TENSOR<T,0,2>& ddstrain,
        const MATRIX<T,2,1>& rdstrain,const MATRIX<T,1>& rxstrain,const VECTOR<T,2>& sigma){
    dPi_dF.Set_Hessian_Block(Contract(ddstrain,P.x)+SYMMETRIC_MATRIX<T,2>::Conjugate_With_Transpose(dstrain,dPi_dF.Get_Hessian_Block()));
    T c(dPi_dF.x1010+dPi_dF.x1001);
    T c_hat=rdstrain.Column(0).Dot(P.x)+rxstrain(0,0)*c;
    T d_hat((P_hat(1)+P_hat(0))/(sigma(1)+sigma(0))); 
    dPi_dF.x1010=0.5*(c_hat+d_hat);
    dPi_dF.x1001=0.5*(c_hat-d_hat);
}
template<class T> void
Transform_Isotropic_Stress_Derivative(DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,
        const DIAGONAL_MATRIX<T,3>& P,const DIAGONAL_MATRIX<T,3>& P_hat,
        const MATRIX<T,3>& dstrain,const SYMMETRIC_TENSOR<T,0,3>& ddstrain,
        const MATRIX<T,3,3>& rdstrain,const MATRIX<T,3>& rxstrain,const VECTOR<T,3>& sigma){
    dPi_dF.Set_Hessian_Block(Contract(ddstrain,P.x)+SYMMETRIC_MATRIX<T,3>::Conjugate_With_Transpose(dstrain,dPi_dF.Get_Hessian_Block()));
    VECTOR<T,3> c(dPi_dF.x2112+dPi_dF.x2121,dPi_dF.x2020+dPi_dF.x2002,dPi_dF.x1010+dPi_dF.x1001);
    VECTOR<T,3> c_hat=rdstrain.Transpose_Times(P.x)+rxstrain.Transpose_Times(c);
    VECTOR<T,3> d_hat((P_hat(1)+P_hat(2))/(sigma(1)+sigma(2)),(P_hat(0)+P_hat(2))/(sigma(0)+sigma(2)),(P_hat(1)+P_hat(0))/(sigma(1)+sigma(0))); 
    dPi_dF.x1010=0.5*(c_hat(2)+d_hat(2));
    dPi_dF.x1001=0.5*(c_hat(2)-d_hat(2));
    dPi_dF.x2020=0.5*(c_hat(1)+d_hat(1));
    dPi_dF.x2002=0.5*(c_hat(1)-d_hat(1));
    dPi_dF.x2121=0.5*(c_hat(0)+d_hat(0));
    dPi_dF.x2112=0.5*(c_hat(0)-d_hat(0));
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

    if(plasticity && plasticity->use_implicit){
#pragma omp parallel for
        for(int k=0;k<gather_scatter.simulated_particles.m;k++){
            int p=gather_scatter.simulated_particles(k);
            MATRIX<T,TV::m> V_local;
            particles.F(p).Fast_Singular_Value_Decomposition(U(p),sigma(p),V_local);
            if(want_dE) FV(p)=force_helper.Fn(p)*V_local;
            TV strain;
            MATRIX<T,TV::m> dstrain;
            SYMMETRIC_TENSOR<T,0,TV::m> ddstrain;
            MATRIX<T,TV::m,TV::SPIN::m> rdstrain;
            MATRIX<T,TV::SPIN::m> rxstrain;
            if(plasticity->Compute(strain,want_dE?&dstrain:0,want_ddE?&ddstrain:0,
                        want_ddE?&rdstrain:0,want_ddE?&rxstrain:0,sigma(p).x,false,p)){
                if(want_dE){
                    DIAGONAL_MATRIX<T,TV::m> P_original=constitutive_model.P_From_Strain(DIAGONAL_MATRIX<T,TV::m>(strain),p);
                    DIAGONAL_MATRIX<T,TV::m> P_hat(dstrain.Transpose_Times(P_original.x));
                    PFT(p)=U(p)*P_hat.Times_Transpose(FV(p));
                    if(want_ddE){
                        constitutive_model.Isotropic_Stress_Derivative(DIAGONAL_MATRIX<T,TV::m>(strain),dPi_dF(p),p);
                        Transform_Isotropic_Stress_Derivative(dPi_dF(p),P_original,P_hat,dstrain,ddstrain,rdstrain,rxstrain,sigma(p).x);
                    }
                }
                sigma(p).x=strain;
            }
            else{
                if(want_dE) PFT(p)=U(p)*constitutive_model.P_From_Strain(sigma(p),p).Times_Transpose(FV(p));
                if(want_ddE) constitutive_model.Isotropic_Stress_Derivative(sigma(p),dPi_dF(p),p);}}}
    else{
#pragma omp parallel for
        for(int k=0;k<gather_scatter.simulated_particles.m;k++){
            int p=gather_scatter.simulated_particles(k);
            MATRIX<T,TV::m> V_local;
            particles.F(p).Fast_Singular_Value_Decomposition(U(p),sigma(p),V_local);
            if(want_dE){
                FV(p)=force_helper.Fn(p)*V_local;
                PFT(p)=U(p)*constitutive_model.P_From_Strain(sigma(p),p).Times_Transpose(FV(p));}
            if(want_ddE) constitutive_model.Isotropic_Stress_Derivative(sigma(p),dPi_dF(p),p);}}
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
