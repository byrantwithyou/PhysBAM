//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/Robust_Functions.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ST_VENANT_KIRCHHOFF_HENCKY_STRAIN.h>
#include <cfloat>
using namespace PhysBAM;
using ::std::log;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
ST_VENANT_KIRCHHOFF_HENCKY_STRAIN(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient)
    :youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input)
{
    assert(-1<poissons_ratio && poissons_ratio<.5);
    constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    constant_mu=youngs_modulus/(2*(1+poissons_ratio));
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;
    //if(d==2) failure_threshold=sqrt((constant_mu+constant_lambda)/(3*constant_mu+3*constant_lambda));
    //else failure_threshold=sqrt((constant_mu+(T)1.5*constant_lambda)/(3*constant_mu+(T)4.5*constant_lambda));
    failure_threshold=0;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
~ST_VENANT_KIRCHHOFF_HENCKY_STRAIN()
{
}
//#####################################################################
// Function Zero_Out_Mu
//#####################################################################
template<class T,int d> void ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
Zero_Out_Mu()
{
    constant_mu=0;
    constant_beta=0;
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    if(F.x(d-1)<=0)
        return FLT_MAX;
    else{
        T id_mu=Mu(id),id_lambda=Lambda(id);
        DIAGONAL_MATRIX<T,d> log_F=log(F);
        T sum_sqr=sqr(log_F.Trace());
        return id_mu*log_F.Frobenius_Norm_Squared()+(T)0.5*id_lambda*sum_sqr;}
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    if(F.x(d-1)<=0)
        return DIAGONAL_MATRIX<T,d>();
    else{
        T id_mu=Mu(id),id_lambda=Lambda(id);
        DIAGONAL_MATRIX<T,d> log_E=log(F);
        DIAGONAL_MATRIX<T,d> F_Inv=F.Inverse();
        return ((T)2*id_mu*log_E+id_lambda*log_E.Trace())*F_Inv;}
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    return DIAGONAL_MATRIX<T,d>();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dPi_dF,T failure_threshold,T mu,T lambda)
{
    if(F.x.y<=0)
        dPi_dF.x0000=dPi_dF.x1111=dPi_dF.x1100=dPi_dF.x1010=dPi_dF.x1001=(T)0;
    else{
        DIAGONAL_MATRIX<T,2> F_threshold=F.Clamp_Min(failure_threshold);
        T lambda_plus_two_mu=lambda+2*mu,two_mu=2*mu,s0=F_threshold.x.x,s1=F_threshold.x.y;
        T l0=log(s0),l1=log(s1),sum_log=l0+l1;
        SYMMETRIC_MATRIX<T,2> F_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(VECTOR<T,2>(F_threshold.x.x,F_threshold.x.y));
        T Psi_0_minus_Psi_1=(two_mu/s0)*(diff_log_over_diff(s0,s1)-l1/s1)-lambda*sum_log/F_outer.x10;
        T Psi_0_plus_Psi_1=((lambda*s0+lambda_plus_two_mu*s1)*l0+(lambda*s1+lambda_plus_two_mu*s0)*l1)/(F_outer.x10*(s0+s1));
        dPi_dF.x0000=(lambda_plus_two_mu*(1-l0)-lambda*l1)/F_outer.x00;
        dPi_dF.x1111=(lambda_plus_two_mu*(1-l1)-lambda*l0)/F_outer.x11;
        dPi_dF.x1100=lambda/F_outer.x10;
        dPi_dF.x1010=(T)0.5*(Psi_0_minus_Psi_1+Psi_0_plus_Psi_1);
        dPi_dF.x1001=(T)0.5*(Psi_0_minus_Psi_1-Psi_0_plus_Psi_1);}
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,T failure_threshold,T mu,T lambda)
{
    if(F.x.z<=0){
        dPi_dF.x0000=dPi_dF.x1111=dPi_dF.x2222=dPi_dF.x1100=dPi_dF.x2200=dPi_dF.x2211=dPi_dF.x1010=dPi_dF.x2020=dPi_dF.x2121=dPi_dF.x1001=dPi_dF.x2002=dPi_dF.x2112=(T)0;}
    else{
        DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold);
        T lambda_plus_two_mu=lambda+2*mu,two_mu=2*mu,s0=F_threshold.x.x,s1=F_threshold.x.y,s2=F_threshold.x.z;
        T l0=log(s0),l1=log(s1),l2=log(s2),sum_log=l0+l1+l2;
        SYMMETRIC_MATRIX<T,3> F_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(VECTOR<T,3>(s0,s1,s2));
        MATRIX<T,3> sl_outer=MATRIX<T,3>::Outer_Product(VECTOR<T,3>(s0,s1,s2),VECTOR<T,3>(l0,l1,l2));
        //Hessian variables
        dPi_dF.x0000=(lambda_plus_two_mu-two_mu*l0-lambda*sum_log)/F_outer.x00;
        dPi_dF.x1111=(lambda_plus_two_mu-two_mu*l1-lambda*sum_log)/F_outer.x11;
        dPi_dF.x2222=(lambda_plus_two_mu-two_mu*l2-lambda*sum_log)/F_outer.x22;
        dPi_dF.x1100=lambda/F_outer.x10;
        dPi_dF.x2200=lambda/F_outer.x20;
        dPi_dF.x2211=lambda/F_outer.x21;
        //2x2 block-matrices
        T Psi_0_minus_Psi_1=(two_mu/s0)*(diff_log_over_diff(s0,s1)-l1/s1)-lambda*sum_log/F_outer.x10;
        T Psi_0_minus_Psi_2=(two_mu/s0)*(diff_log_over_diff(s0,s2)-l2/s2)-lambda*sum_log/F_outer.x20;
        T Psi_1_minus_Psi_2=(two_mu/s1)*(diff_log_over_diff(s1,s2)-l2/s2)-lambda*sum_log/F_outer.x21;
        T Psi_0_plus_Psi_1=(two_mu*(s1*l0+s0*l1)+lambda*(sl_outer.Row(0).Sum()+sl_outer.Row(1).Sum()))/(F_outer.x10*(s0+s1));
        T Psi_0_plus_Psi_2=(two_mu*(s2*l0+s0*l2)+lambda*(sl_outer.Row(0).Sum()+sl_outer.Row(2).Sum()))/(F_outer.x20*(s0+s2));
        T Psi_1_plus_Psi_2=(two_mu*(s2*l1+s1*l2)+lambda*(sl_outer.Row(1).Sum()+sl_outer.Row(2).Sum()))/(F_outer.x21*(s1+s2));
        dPi_dF.x1010=(T)0.5*(Psi_0_minus_Psi_1+Psi_0_plus_Psi_1);
        dPi_dF.x1001=(T)0.5*(Psi_0_minus_Psi_1-Psi_0_plus_Psi_1); 
        dPi_dF.x2020=(T)0.5*(Psi_0_minus_Psi_2+Psi_0_plus_Psi_2);
        dPi_dF.x2002=(T)0.5*(Psi_0_minus_Psi_2-Psi_0_plus_Psi_2); 
        dPi_dF.x2121=(T)0.5*(Psi_1_minus_Psi_2+Psi_1_plus_Psi_2);
        dPi_dF.x2112=(T)0.5*(Psi_1_minus_Psi_2-Psi_1_plus_Psi_2);}
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    Isotropic_Stress_Derivative_Helper(F,dP_dF,failure_threshold,id_mu,id_lambda);
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
namespace PhysBAM{
template class ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<double,2>;
template class ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<double,3>;
template class ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<float,2>;
template class ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<float,3>;
}
