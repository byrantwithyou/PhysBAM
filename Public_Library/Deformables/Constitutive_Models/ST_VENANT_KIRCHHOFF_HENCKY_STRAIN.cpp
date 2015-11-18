//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ST_VENANT_KIRCHHOFF_HENCKY_STRAIN.h>
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
    if(d==2) failure_threshold=sqrt((constant_mu+constant_lambda)/(3*constant_mu+3*constant_lambda));
    else failure_threshold=sqrt((constant_mu+(T)1.5*constant_lambda)/(3*constant_mu+(T)4.5*constant_lambda));
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
    T id_mu=(mu.m?mu(id):constant_mu),id_lambda=(lambda.m?lambda(id):constant_lambda);
    DIAGONAL_MATRIX<T,d> log_F=log(F.Abs());
    T sum_sqr=sqr(log_F.Trace());
    return id_mu*log_F.Frobenius_Norm_Squared()+(T)0.5*id_lambda*sum_sqr;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=(mu.m?mu(id):constant_mu),id_lambda=(lambda.m?lambda(id):constant_lambda);
    DIAGONAL_MATRIX<T,d> log_E=log(F.Abs());
    DIAGONAL_MATRIX<T,d> F_Inv=F.Inverse();
    return ((T)2*id_mu*log_E+id_lambda*log_E.Trace())*F_Inv;
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return DIAGONAL_MATRIX<T,d>();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dPi_dF,T failure_threshold,T mu,T lambda)
{
    DIAGONAL_MATRIX<T,2> F_threshold=F.Clamp_Min(failure_threshold);
    T lambda_plus_two_mu=lambda+2*mu,two_mu=2*mu,two_lambda_plus_two_mu=2*(lambda+mu),s0=F_threshold.x.x,s1=F_threshold.x.y;
    T l0=log(s0),l1=log(s1);
    T s0_sqr=sqr(s0),s1_sqr=sqr(s1);
    T eps=s1-s0;
    dPi_dF.x0000=(lambda_plus_two_mu*(1-l0)-lambda*l1)/s0_sqr;
    dPi_dF.x1111=(lambda_plus_two_mu*(1-l1)-lambda*l0)/s1_sqr;
    dPi_dF.x1100=lambda/(s0*s1);
    dPi_dF.x1010=(abs(eps)<1e-8)?mu/s0_sqr:two_mu*(l0-l1)/(s0_sqr-s1_sqr);
    dPi_dF.x1001=(abs(eps)<1e-8)?(mu-two_lambda_plus_two_mu*l0)/s0_sqr:((-lambda*s0_sqr+lambda_plus_two_mu*s1_sqr)*l0+(lambda*s1_sqr-lambda_plus_two_mu*s0_sqr)*l1)/(s0*s1*(s0_sqr-s1_sqr));
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,T failure_threshold,T mu,T lambda)
{
    DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold);
    T lambda_plus_two_mu=lambda+2*mu,two_mu=2*mu,two_lambda_plus_two_mu=2*(lambda+mu),s0=F_threshold.x.x,s1=F_threshold.x.y,s2=F_threshold.x.y;
    T s0_sqr=sqr(s0),s1_sqr=sqr(s1),s2_sqr=sqr(s2);
    T eps_10=s1-s0,eps_20=s2-s0,eps_21=s2-s1;
    T l0=log(s0),l1=log(s1),l2=log(s2),lt=log(s0*s1*s2);
    //Hessian variables
    dPi_dF.x0000=(lambda_plus_two_mu-two_mu*l0-lambda*lt)/s0_sqr;
    dPi_dF.x1111=(lambda_plus_two_mu-two_mu*l1-lambda*lt)/s1_sqr;
    dPi_dF.x2222=(lambda_plus_two_mu-two_mu*l2-lambda*lt)/s2_sqr;
    dPi_dF.x1100=lambda/(s0*s1);
    dPi_dF.x2200=lambda/(s0*s2);
    dPi_dF.x2211=lambda/(s1*s2);
    //ijij
    dPi_dF.x1010=abs(eps_10<1e-8)?mu/s0_sqr:two_mu*(l0-l1)/(s0_sqr-s1_sqr);
    dPi_dF.x2020=abs(eps_20<1e-8)?mu/s0_sqr:two_mu*(l0-l2)/(s0_sqr-s2_sqr);
    dPi_dF.x2121=abs(eps_21<1e-8)?mu/s1_sqr:two_mu*(l1-l2)/(s1_sqr-s2_sqr);
    //ijji
    dPi_dF.x1001=abs(eps_10<1e-8)?(mu-two_lambda_plus_two_mu*l0-lambda*l2)/s0_sqr:((-s0_sqr*lambda+s1_sqr*lambda_plus_two_mu)*l0+(s1_sqr*lambda-s0_sqr*lambda_plus_two_mu)*l1+(-s0_sqr+s1_sqr)*lambda*l2)/(s0*s1*(s0_sqr-s1_sqr));
    dPi_dF.x2002=abs(eps_20<1e-8)?(mu-two_lambda_plus_two_mu*l0-lambda*l1)/s0_sqr:((-s0_sqr*lambda+s2_sqr*lambda_plus_two_mu)*l0+(-s0_sqr+s2_sqr)*lambda*l1+(s2_sqr*lambda-s0_sqr*lambda_plus_two_mu)*l2)/(s0*s2*(s0_sqr-s2_sqr));
    dPi_dF.x2112=abs(eps_21<1e-8)?(mu-two_lambda_plus_two_mu*l1-lambda*l0)/s1_sqr:((-s1_sqr+s2_sqr)*lambda*l0+(-s1_sqr*lambda+s2_sqr*lambda_plus_two_mu)*l1+(s2_sqr*lambda-s1_sqr*lambda_plus_two_mu)*l2)/(s1*s2*(s1_sqr-s2_sqr));
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int id) const
{
    T id_mu=(mu.m?mu(id):constant_mu),id_lambda=(lambda.m?lambda(id):constant_lambda);
    Isotropic_Stress_Derivative_Helper(F,dP_dF,failure_threshold,id_mu,id_lambda);
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
namespace PhysBAM{
template class ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<double,2>;
template class ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<double,3>;
template class ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<float,2>;
template class ST_VENANT_KIRCHHOFF_HENCKY_STRAIN<float,3>;
}
