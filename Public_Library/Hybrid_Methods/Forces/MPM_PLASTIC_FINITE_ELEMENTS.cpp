//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/MPM_PLASTIC_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_PLASTICITY_MODEL.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
namespace PhysBAM{
//#####################################################################
// Function Times
//#####################################################################
template<class TV> auto MPM_PLASTIC_FINITE_ELEMENTS<TV>::HESSIAN_HELPER::
Times(const MATRIX<T,TV::m>& M) const -> MATRIX<T,TV::m>
{
    MATRIX<T,TV::m> A(DIAGONAL_MATRIX<T,TV::m>(diag*M.Diagonal_Part().x));
    for(int k=0;k<TV::SPIN::m;k++){
        int i=(k+1)%TV::m;
        int j=(k+2)%TV::m;
        T mij=M(i,j),mji=M(j,i);
        T d=off_diag_diff(k)*(mij+mji);
        T s=off_diag_sum(k)*(mij-mji);
        A(i,j)=(d+s)/2;
        A(j,i)=(d-s)/2;}
    return A;
}
//#####################################################################
// Function Set
//#####################################################################
template<class TV> void MPM_PLASTIC_FINITE_ELEMENTS<TV>::HESSIAN_HELPER::
Set(const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dPi_dF)
{
    diag=dPi_dF.H;
    off_diag_diff=dPi_dF.B+dPi_dF.C;
    off_diag_sum=dPi_dF.B-dPi_dF.C;
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PLASTIC_FINITE_ELEMENTS<TV>::
MPM_PLASTIC_FINITE_ELEMENTS(const MPM_FORCE_HELPER<TV>& force_helper,
    ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m>& constitutive_model,
    GATHER_SCATTER<TV>& gather_scatter_input,ARRAY<int>* affected_particles,
    MPM_PLASTICITY_MODEL<TV>& plasticity)
    :BASE(force_helper),constitutive_model(constitutive_model),plasticity(plasticity),affect_all(!affected_particles),
    gather_scatter(affect_all?gather_scatter_input:*new GATHER_SCATTER<TV>(gather_scatter_input.grid,*new ARRAY<int>(*affected_particles)))
{
    if(!affect_all){
        gather_scatter.weights=gather_scatter_input.weights;
        gather_scatter.threads=gather_scatter_input.threads;
        gather_scatter.Prepare_Scatter(particles);}
    plasticity.gather_scatter=&gather_scatter;
    constitutive_model.mu=&particles.mu;
    constitutive_model.lambda=&particles.lambda;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PLASTIC_FINITE_ELEMENTS<TV>:: 
~MPM_PLASTIC_FINITE_ELEMENTS()
{
    if(!affect_all){
        delete &gather_scatter.simulated_particles;
        delete &gather_scatter;}
    delete &constitutive_model;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void MPM_PLASTIC_FINITE_ELEMENTS<TV>:: 
Precompute(const T time,const T dt,bool want_dE,bool want_ddE)
{
    PHYSBAM_ASSERT(want_dE);
    U.Resize(particles.X.m);
    FV.Resize(particles.X.m);
    hessian_helper.Resize(particles.X.m);
    PFT.Resize(particles.X.m);

    int num_projected=0;
    if(plasticity.use_implicit){
#pragma omp parallel for reduction(+:num_projected)
        for(int k=0;k<gather_scatter.simulated_particles.m;k++){
            int p=gather_scatter.simulated_particles(k);
            MATRIX<T,TV::m> V_local;
            DIAGONAL_MATRIX<T,TV::m> sigma;
            particles.F(p).Singular_Value_Decomposition(U(p),sigma,V_local);
            FV(p)=force_helper.Fn(p)*V_local;
            TV strain;
            MATRIX<T,TV::m> dstrain;
            typename TV::SPIN r_sum,r_diff;
            MATRIX<T,TV::m,TV::SPIN::m> rdstrain;
            MATRIX<T,TV::SPIN::m> rxstrain;
            bool projected=plasticity.Compute(strain,want_ddE?&dstrain:0,want_ddE?&r_sum:0,want_ddE?&r_diff:0,sigma.x,false,p);
            num_projected+=projected;
            DIAGONAL_MATRIX<T,TV::m> new_sigma=projected?DIAGONAL_MATRIX<T,TV::m>(strain):sigma;
            DIAGONAL_MATRIX<T,TV::m> P_hat=constitutive_model.P_From_Strain(new_sigma,p);
            PFT(p)=U(p)*P_hat.Times_Transpose(FV(p));

            if(want_ddE){
                DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV> dPi_dF;
                constitutive_model.Isotropic_Stress_Derivative(new_sigma,dPi_dF,p);
                HESSIAN_HELPER& hh=hessian_helper(p);
                hh.Set(dPi_dF);
                if(projected){
                    hh.diag=hh.diag*dstrain;
                    hh.off_diag_diff=hh.off_diag_diff*r_diff;
                    hh.off_diag_sum=hh.off_diag_sum*r_sum;}}}}
    LOG::printf("PROJECTED: %P\n",num_projected);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR MPM_PLASTIC_FINITE_ELEMENTS<TV>:: 
Potential_Energy(const T time) const
{
    return 0;
}
//#####################################################################
// Function Add_Forces
//#####################################################################
template<class TV> void MPM_PLASTIC_FINITE_ELEMENTS<TV>:: 
Add_Forces(ARRAY<TV,TV_INT>& F,const T time) const
{
    T c=force_helper.quad_F_coeff;
    bool use_c=c!=0;
    gather_scatter.template Scatter<MATRIX<T,TV::m> >(true,
        [this,c,use_c](int p,MATRIX<T,TV::m>& AF)
        {
            MATRIX<T,TV::m> P=PFT(p);
            if(use_c){MATRIX<T,TV::m> B=force_helper.B(p);P+=(P.Times_Transpose(B)+B.Transpose_Times(P))*c;}
            AF=P*particles.volume(p);
        },
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,const MATRIX<T,TV::m>& A)
        {F(it.Index())-=A*it.Gradient();});
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_PLASTIC_FINITE_ELEMENTS<TV>:: 
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& Z,const T time) const
{
    tmp.Resize(particles.X.m);
    
    T c=force_helper.quad_F_coeff;
    bool use_c=c!=0;
    gather_scatter.template Gather<int>(true,
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
            MATRIX<T,TV::m> M=UU*hessian_helper(p).Times(G).Times_Transpose(FVV);
            if(use_c){
                MATRIX<T,TV::m> B=force_helper.B(p),P=PFT(p);
                M+=(M.Times_Transpose(B)+B.Transpose_Times(M)+P.Times_Transpose(tmp(p))+tmp(p).Transpose_Times(P))*c;}
            tmp(p)=M*particles.volume(p);
        });

    gather_scatter.template Scatter<int>(true,
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {F(it.Index())+=tmp(p)*it.Gradient();});
}
template class MPM_PLASTIC_FINITE_ELEMENTS<VECTOR<float,2> >;
template class MPM_PLASTIC_FINITE_ELEMENTS<VECTOR<float,3> >;
template class MPM_PLASTIC_FINITE_ELEMENTS<VECTOR<double,2> >;
template class MPM_PLASTIC_FINITE_ELEMENTS<VECTOR<double,3> >;
}
