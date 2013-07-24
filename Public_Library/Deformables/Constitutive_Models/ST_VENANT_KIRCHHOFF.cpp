//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Russell Howes, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ST_VENANT_KIRCHHOFF.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> ST_VENANT_KIRCHHOFF<T,d>::
ST_VENANT_KIRCHHOFF(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient)
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
template<class T,int d> ST_VENANT_KIRCHHOFF<T,d>::
~ST_VENANT_KIRCHHOFF()
{
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T ST_VENANT_KIRCHHOFF<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    DIAGONAL_MATRIX<T,d> F_threshold=F.Clamp_Min(failure_threshold),strain=(F_threshold*F_threshold-1)*(T).5,strain_squared=strain*strain;
    return (T).5*constant_lambda*strain.Trace()*strain.Trace()+constant_mu*strain_squared.Trace();
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> ST_VENANT_KIRCHHOFF<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    DIAGONAL_MATRIX<T,d> F_threshold=F.Clamp_Min(failure_threshold),twice_strain=F_threshold*F_threshold-1;
    return F_threshold*(scale*constant_mu*twice_strain+(T).5*scale*constant_lambda*twice_strain.Trace());
}
template<class T,int d> MATRIX<T,d> ST_VENANT_KIRCHHOFF<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    DIAGONAL_MATRIX<T,d> F_threshold=F.Clamp_Min(failure_threshold);
    SYMMETRIC_MATRIX<T,d> strain_rate=(F_threshold*F_dot).Symmetric_Part();
    return F_threshold*(2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace());
}
template<class T> static void Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,
    T failure_threshold,T mu,T lambda)
{
    DIAGONAL_MATRIX<T,2> F_threshold=F.Clamp_Min(failure_threshold);
    T lambda_tr_G_minus_mu=(T).5*lambda*(F_threshold*F_threshold-1).Trace()-mu,three_mu_plus_lambda=3*mu+lambda;
    SYMMETRIC_MATRIX<T,2> F_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(VECTOR<T,2>(F_threshold.x11,F_threshold.x22));
    dP_dF.x1111=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x11;//alpha+beta+gamma
    dP_dF.x2222=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x22;
    dP_dF.x2211=lambda*F_outer.x21;//gamma
    dP_dF.x2121=lambda_tr_G_minus_mu+mu*(F_outer.x22+F_outer.x11);//alpha
    dP_dF.x2112=mu*F_outer.x21;//beta
}
template<class T> static void Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,
    T failure_threshold,T mu,T lambda)
{
    DIAGONAL_MATRIX<T,3> F_threshold=F.Clamp_Min(failure_threshold);
    T lambda_tr_G_minus_mu=(T).5*lambda*(F_threshold*F_threshold-1).Trace()-mu,three_mu_plus_lambda=3*mu+lambda;
    SYMMETRIC_MATRIX<T,3> F_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(VECTOR<T,3>(F_threshold.x11,F_threshold.x22,F_threshold.x33));
    //alpha+beta+gamma
    dPi_dF.x1111=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x11;
    dPi_dF.x2222=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x22;
    dPi_dF.x3333=lambda_tr_G_minus_mu+three_mu_plus_lambda*F_outer.x33;
    //gamma
    dPi_dF.x2211=lambda*F_outer.x21;
    dPi_dF.x3311=lambda*F_outer.x31;
    dPi_dF.x3322=lambda*F_outer.x32;
    //alpha
    dPi_dF.x2121=lambda_tr_G_minus_mu+mu*(F_outer.x22+F_outer.x11);
    dPi_dF.x3131=lambda_tr_G_minus_mu+mu*(F_outer.x33+F_outer.x11);
    dPi_dF.x3232=lambda_tr_G_minus_mu+mu*(F_outer.x33+F_outer.x22);
    //beta
    dPi_dF.x2112=mu*F_outer.x21;
    dPi_dF.x3113=mu*F_outer.x31;
    dPi_dF.x3223=mu*F_outer.x32;
}
template<class T,int d> void ST_VENANT_KIRCHHOFF<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int simplex) const
{
    Isotropic_Stress_Derivative_Helper(F,dP_dF,failure_threshold,constant_mu,constant_lambda);
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
namespace PhysBAM{
template class ST_VENANT_KIRCHHOFF<double,2>;
template class ST_VENANT_KIRCHHOFF<double,3>;
template class ST_VENANT_KIRCHHOFF<float,2>;
template class ST_VENANT_KIRCHHOFF<float,3>;
}
