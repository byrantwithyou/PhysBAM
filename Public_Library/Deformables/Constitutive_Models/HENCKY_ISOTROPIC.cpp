//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Russell Howes, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/HENCKY_ISOTROPIC.h>
using namespace PhysBAM;
using ::std::log;
//#####################################################################
// Constructor
// TODO: Precompute log F?
//#####################################################################
template<class T,int d> HENCKY_ISOTROPIC<T,d>::
HENCKY_ISOTROPIC(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient)
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
template<class T,int d> HENCKY_ISOTROPIC<T,d>::
~HENCKY_ISOTROPIC()
{
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T HENCKY_ISOTROPIC<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=(mu.m?mu(id):constant_mu),id_lambda=(lambda.m?lambda(id):constant_lambda);
    DIAGONAL_MATRIX<T,d> log_F=log(F);
    T tr_E=log_F.Trace();
    T E_sq=(log_F*log_F).Trace();
    return id_mu*E_sq+(T)0.5*id_lambda*tr_E*tr_E;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> HENCKY_ISOTROPIC<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=(mu.m?mu(id):constant_mu),id_lambda=(lambda.m?lambda(id):constant_lambda);
    DIAGONAL_MATRIX<T,d> log_E=log(F);
    DIAGONAL_MATRIX<T,d> F_Inv=F.Inverse();
    return 2*id_mu*F_Inv*log_E+id_lambda*log_E.Trace()*F_Inv;
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> HENCKY_ISOTROPIC<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return 0;
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,
    T failure_threshold,T mu,T lambda)
{
    DIAGONAL_MATRIX<T,2> F_threshold=F.Clamp_Min(failure_threshold);
    T lambda_tr_G_minus_mu=(T).5*lambda*(F_threshold*F_threshold-1).Trace()-mu,three_mu_plus_lambda=3*mu+lambda;
    SYMMETRIC_MATRIX<T,2> F_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(VECTOR<T,2>(F_threshold.x.x,F_threshold.x.y));
    dP_dF.x0000=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x00;//alpha+beta+gamma
    dP_dF.x1111=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x11;
    dP_dF.x1100=lambda*F_outer.x10;//gamma
    dP_dF.x1010=lambda_tr_G_minus_mu+mu*(F_outer.x11+F_outer.x00);//alpha
    dP_dF.x1001=mu*F_outer.x10;//beta
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,
    T failure_threshold,T mu,T lambda)
{
    DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold);
    T lambda_tr_G_minus_mu=(T).5*lambda*(F_threshold*F_threshold-1).Trace()-mu,three_mu_plus_lambda=3*mu+lambda;
    SYMMETRIC_MATRIX<T,3> F_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(VECTOR<T,3>(F_threshold.x.x,F_threshold.x.y,F_threshold.x.z));
    //alpha+beta+gamma
    dPi_dF.x0000=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x00;
    dPi_dF.x1111=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x11;
    dPi_dF.x2222=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x22;
    //gamma
    dPi_dF.x1100=lambda*F_outer.x10;
    dPi_dF.x2200=lambda*F_outer.x20;
    dPi_dF.x2211=lambda*F_outer.x21;
    //alpha
    dPi_dF.x1010=lambda_tr_G_minus_mu+mu*(F_outer.x11+F_outer.x00);
    dPi_dF.x2020=lambda_tr_G_minus_mu+mu*(F_outer.x22+F_outer.x00);
    dPi_dF.x2121=lambda_tr_G_minus_mu+mu*(F_outer.x22+F_outer.x11);
    //beta
    dPi_dF.x1001=mu*F_outer.x10;
    dPi_dF.x2002=mu*F_outer.x20;
    dPi_dF.x2112=mu*F_outer.x21;
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void HENCKY_ISOTROPIC<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int id) const
{
    T id_mu=(mu.m?mu(id):constant_mu),id_lambda=(lambda.m?lambda(id):constant_lambda);


    Isotropic_Stress_Derivative_Helper(F,dP_dF,failure_threshold,id_mu,id_lambda);
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
namespace PhysBAM{
template class HENCKY_ISOTROPIC<double,2>;
template class HENCKY_ISOTROPIC<double,3>;
template class HENCKY_ISOTROPIC<float,2>;
template class HENCKY_ISOTROPIC<float,3>;
}
