//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/Robust_Functions.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Constitutive_Models/MPM_DRUCKER_PRAGER.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_DRUCKER_PRAGER<TV>::
MPM_DRUCKER_PRAGER(MPM_PARTICLES<TV>& particles,T a0,T a1,T a3,T a4)
    :MPM_PLASTICITY_MODEL<TV>(particles),a0(a0),a1(a1),a3(a3),a4(a4)
{
    particles.Add_Array(ATTRIBUTE_ID_PLASTIC_DEFORMATION,&plastic_def);
    particles.Add_Array(ATTRIBUTE_ID_DP_RHO_F,&rho_F);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_DRUCKER_PRAGER<TV>::
~MPM_DRUCKER_PRAGER()
{
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> void MPM_DRUCKER_PRAGER<TV>::
Initialize_Particle(int p) const
{
    Update_Hardening(p,0);
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> bool MPM_DRUCKER_PRAGER<TV>::
Compute(TV& strain,MATRIX<T,TV::m>* dstrain,SYMMETRIC_TENSOR<T,0,TV::m>* ddstrain,
    MATRIX<T,TV::m,TV::SPIN::m>* rdstrain,MATRIX<T,TV::SPIN::m>* rxstrain,const TV& Fe,bool store_hardening,int id) const
{
    typedef typename TV::SPIN TVSPIN;
    T mu=particles.mu(id),lambda=particles.lambda(id);
    T g=1/(d*lambda+2*mu);
    T b=1/(2*mu);
    TV strain_trial=log(abs(Fe));
    T k=strain_trial.Sum();
    TV sh=strain_trial-k/d;
    T q=sh.Magnitude();
    T r=rho_F(id)*k/g;
    T dg=q+b*r;
    if(q+g*g*r/b<=0) return false;
    if(q==0 || k>0){
        strain=TV::All_Ones_Vector();
        if(dstrain) *dstrain=MATRIX<T,d>();
        if(ddstrain) *ddstrain=SYMMETRIC_TENSOR<T,0,d>();
        if(store_hardening) Update_Hardening(id,strain_trial.Magnitude());
        return true;} 
    T p=r*b/q;
    TV h_strain=k/d-p*sh;
    strain=exp(h_strain);
    if(store_hardening) Update_Hardening(id,dg);

// What if Fe<0?
    if(!dstrain && !ddstrain && !rdstrain) return true;
    TV o=TV::All_Ones_Vector();
    DIAGONAL_MATRIX<T,d> DMS(strain);
    TV one_over_Fe=(T)1/Fe;

    DIAGONAL_MATRIX<T,d> D_strain_trial(one_over_Fe);
    TV D_k=one_over_Fe;
    MATRIX<T,d> H1=MATRIX<T,d>::Outer_Product(o,D_k/d);
    TV D_q=D_strain_trial*sh/q;
    TV D_r=rho_F(id)/g*D_k;
    TV D_p=(b*D_r-p*D_q)/q;
    MATRIX<T,d> D_h_strain=(1+p)*H1-MATRIX<T,d>::Outer_Product(sh,D_p)-p*D_strain_trial;
    if(dstrain) *dstrain=DMS*D_h_strain;

    if(!ddstrain) return true;
    TV D_one_over_Fe=-one_over_Fe*one_over_Fe;
    DIAGONAL_TENSOR<T,d> DD_strain_trial(D_one_over_Fe);
    DIAGONAL_MATRIX<T,d> DD_k(D_one_over_Fe);
    SYMMETRIC_TENSOR<T,0,d> D_H1=Tensor_Product<0>(DD_k/d,o);
    SYMMETRIC_TENSOR<T,0,d> DD_sh=DD_strain_trial-D_H1;
    SYMMETRIC_MATRIX<T,d> DD_q=Contract<2>(DD_strain_trial,sh/q)-DD_k/q-SYMMETRIC_MATRIX<T,d>::Outer_Product(D_k)/(d*q)-SYMMETRIC_MATRIX<T,d>::Outer_Product(D_q)/q;
    DIAGONAL_MATRIX<T,d> DD_r=rho_F(id)/g*DD_k;
    SYMMETRIC_MATRIX<T,d> DD_p=(b*DD_r-SYMMETRIC_MATRIX<T,d>::Symmetric_Outer_Product(D_q,D_p)-p*DD_q)/q;
    *ddstrain=Tensor_Product<0>(DD_k,strain/d)-Tensor_Product<0>(DD_p,DMS*sh)-Contract<0,0>(DD_sh,p*DMS)-Symmetric_Tensor_Product<1,2>(DMS*(D_strain_trial-H1),D_p)
        +Symmetric_Double_Contract_12(DIAGONAL_TENSOR<T,d>(strain/2),D_h_strain,D_h_strain);
    if(rdstrain){
        TV R_Dh=-p*strain*one_over_Fe;
        TVSPIN R_Fi_lnFe,R_k,RX_exp,RX_ln,R_Da,R_Db,R_D;
        for(int i=0;i<TVSPIN::m;i++){
            T x=Fe((i+1)%d),y=Fe((i+2)%d),lnx=strain_trial((i+1)%d);
            R_Fi_lnFe(i)=(diff_log_over_diff(x,y)-lnx/x)/y;
            R_k(i)=-1/(x*y);
            T e1=h_strain((i+1)%d),e2=h_strain((i+2)%d);
            RX_exp(i)=diff_exp_over_diff(e1,e2);
            RX_ln(i)=diff_log_over_diff(x,y);
            R_Da(i)=strain((i+2)%d);
            R_Db(i)=one_over_Fe((i+1)%d);
            R_D(i)=R_Dh((i+1)%d);}
        TVSPIN RX_h_strain=-p*RX_ln;
        TVSPIN RX_strain=RX_exp*RX_h_strain;
        TVSPIN RX_D=-p*(DIAGONAL_MATRIX<T,TVSPIN::m>(R_Da)*R_k+DIAGONAL_MATRIX<T,TVSPIN::m>(R_Db)*RX_strain);
        MATRIX<T,d,TVSPIN::m> RH1=MATRIX<T,d,TVSPIN::m>::Outer_Product(o,R_k/d);
        TVSPIN R_q=R_Fi_lnFe/q-R_k*k/(d*q);
        TVSPIN R_r=rho_F(id)/g*R_k;
        TVSPIN R_p=(b*R_r-p*R_q)/q;
        MATRIX<T,d,TVSPIN::m> R_h_strain=(1+p)*RH1-MATRIX<T,d,TVSPIN::m>::Outer_Product(sh,R_p);
        MATRIX<T,d,TVSPIN::m> A=DMS*R_h_strain; // + (D*R-R*E)
        for(int i=0;i<TVSPIN::m;i++)
            A((i+2)%d,i%d)+=RX_D(i);
        *rdstrain=A;
        *rxstrain=DIAGONAL_MATRIX<T,TVSPIN::m>(R_D);
    }

    return true;
}
//#####################################################################
// Function Project_Stress
//#####################################################################
template<class TV> void MPM_DRUCKER_PRAGER<TV>::
Update_Hardening(int id,T plastic_def_increment) const
{
    plastic_def(id)+=plastic_def_increment;
    T phi_F=pi/180*(a0+(a1*plastic_def(id)-a4)*exp(-a3*plastic_def(id)));
    T sin_phi_F=sin(phi_F);
    rho_F(id)=sqrt(2.0/3.0)*2*sin_phi_F/(3-sin_phi_F);
}

static int Register_Attributes()
{
    Register_Attribute_Name(ATTRIBUTE_ID_PLASTIC_DEFORMATION,"plastic_def");
    Register_Attribute_Name(ATTRIBUTE_ID_DP_RHO_F,"rho_F");
    return 0;
}
int MPM_DRUCKER_PRAGER_zzz=Register_Attributes();
template class MPM_DRUCKER_PRAGER<VECTOR<float,2>>;
template class MPM_DRUCKER_PRAGER<VECTOR<float,3>>;
template class MPM_DRUCKER_PRAGER<VECTOR<double,2>>;
template class MPM_DRUCKER_PRAGER<VECTOR<double,3>>;
}
