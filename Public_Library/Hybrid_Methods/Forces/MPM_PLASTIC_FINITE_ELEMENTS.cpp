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
#include <Hybrid_Methods/Forces/MPM_FORCE_HELPER.h>
#include <Hybrid_Methods/Forces/MPM_PLASTIC_FINITE_ELEMENTS.h>
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
// Function Set_Helper
//#####################################################################
template<class T> void
Set_Helper(VECTOR<T,0>& sum,VECTOR<T,0>& diff,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,1>& dPi_dF)
{
}
//#####################################################################
// Function Set_Helper
//#####################################################################
template<class T> void
Set_Helper(VECTOR<T,1>& sum,VECTOR<T,1>& diff,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dPi_dF)
{
    diff(0)=dPi_dF.x1010+dPi_dF.x1001;
    sum(0)=dPi_dF.x1010-dPi_dF.x1001;
}
//#####################################################################
// Function Set_Helper
//#####################################################################
template<class T> void
Set_Helper(VECTOR<T,3>& sum,VECTOR<T,3>& diff,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF)
{
    diff(2)=dPi_dF.x1010+dPi_dF.x1001;
    diff(1)=dPi_dF.x2020+dPi_dF.x2002;
    diff(0)=dPi_dF.x2121+dPi_dF.x2112;
    sum(2)=dPi_dF.x1010-dPi_dF.x1001;
    sum(1)=dPi_dF.x2020-dPi_dF.x2002;
    sum(0)=dPi_dF.x2121-dPi_dF.x2112;
}
//#####################################################################
// Function Set
//#####################################################################
template<class TV> void MPM_PLASTIC_FINITE_ELEMENTS<TV>::HESSIAN_HELPER::
Set(const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,TV::m>& dPi_dF)
{
    diag=dPi_dF.Get_Hessian_Block();
    Set_Helper(off_diag_sum,off_diag_diff,dPi_dF);
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
            particles.F(p).Fast_Singular_Value_Decomposition(U(p),sigma,V_local);
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
                DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,TV::m> dPi_dF;
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
    gather_scatter.template Scatter<MATRIX<T,TV::m> >(
        [this](int p,MATRIX<T,TV::m>& AF)
        {
            MATRIX<T,TV::m> P=PFT(p);
            AF=P*particles.volume(p);
        },
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,const MATRIX<T,TV::m>& A)
        {F(it.Index())-=A*it.Gradient();},
        [](int p,const MATRIX<T,TV::m>& A){},true);
}
//#####################################################################
// Function Add_Hessian_Times
//#####################################################################
template<class TV> void MPM_PLASTIC_FINITE_ELEMENTS<TV>:: 
Add_Hessian_Times(ARRAY<TV,TV_INT>& F,const ARRAY<TV,TV_INT>& Z,const T time) const
{
    tmp.Resize(particles.X.m);
    
    gather_scatter.template Gather<int>(
        [this](int p,int data){tmp(p)=MATRIX<T,TV::m>();},
        [this,&Z](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {
            tmp(p)+=MATRIX<T,TV::m>::Outer_Product(Z(it.Index()),it.Gradient());
        },
        [this](int p,int data){
            MATRIX<T,TV::m> W=tmp(p);
            MATRIX<T,TV::m> UU=U(p),FVV=FV(p),G=UU.Transpose_Times(W*FVV);
            MATRIX<T,TV::m> M=UU*hessian_helper(p).Times(G).Times_Transpose(FVV);
            tmp(p)=M*particles.volume(p);
        },true);

    gather_scatter.template Scatter<int>(
        [this,&F](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
        {F(it.Index())+=tmp(p)*it.Gradient();},true);
}
template class MPM_PLASTIC_FINITE_ELEMENTS<VECTOR<float,2> >;
template class MPM_PLASTIC_FINITE_ELEMENTS<VECTOR<float,3> >;
template class MPM_PLASTIC_FINITE_ELEMENTS<VECTOR<double,2> >;
template class MPM_PLASTIC_FINITE_ELEMENTS<VECTOR<double,3> >;
}
