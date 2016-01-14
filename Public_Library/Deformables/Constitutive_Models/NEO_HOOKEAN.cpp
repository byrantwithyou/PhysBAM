//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/cube.h>
#include <Tools/Math_Tools/pow.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX_1X1.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_1X1.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN<T,d>::
NEO_HOOKEAN(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T failure_threshold_input)
    :youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),failure_threshold(failure_threshold_input),use_constant_ife(false)
{
    assert(poissons_ratio>-1&&poissons_ratio<.5);
    constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    constant_mu=youngs_modulus/(2*(1+poissons_ratio));
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;
    dth_root_failure_threshold=pow<1,d>(failure_threshold);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN<T,d>::
~NEO_HOOKEAN()
{
}
//#####################################################################
// Function Clamp_To_Hyperbola
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,1> NEO_HOOKEAN<T,d>::
Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,1>& F) const
{
    if(sqr(F.x.x)>failure_threshold) return DIAGONAL_MATRIX<T,1>(failure_threshold);
    else return DIAGONAL_MATRIX<T,1>(dth_root_failure_threshold);
}
//#####################################################################
// Function Clamp_To_Hyperbola
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,2> NEO_HOOKEAN<T,d>::
Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,2>& F) const
{
    if(sqr(F.x.x)>failure_threshold) return DIAGONAL_MATRIX<T,2>(F.x.x,failure_threshold/F.x.x);
    else return DIAGONAL_MATRIX<T,2>(dth_root_failure_threshold,dth_root_failure_threshold);
}
//#####################################################################
// Function Clamp_To_Hyperbola
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> NEO_HOOKEAN<T,d>::
Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,3>& F) const
{
    if(F.x.x*F.x.y*F.x.y>failure_threshold) return DIAGONAL_MATRIX<T,3>(F.x.x,F.x.y,failure_threshold/(F.x.x*F.x.y));
    else if(cube(F.x.x)>failure_threshold){
        T clamped=sqrt(failure_threshold/F.x.x);
        return DIAGONAL_MATRIX<T,3>(F.x.x,clamped,clamped);}
    else return DIAGONAL_MATRIX<T,3>(dth_root_failure_threshold,dth_root_failure_threshold,dth_root_failure_threshold);
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
// clamp to hyperbola to avoid indefiniteness "automatically"
template<class T,int d> DIAGONAL_MATRIX<T,d> NEO_HOOKEAN<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id),J=F.Determinant();
    if(J>=failure_threshold) return id_mu*F-(id_mu-id_lambda*log(J))*F.Inverse();
    DIAGONAL_MATRIX<T,d> F_clamp=Clamp_To_Hyperbola(F),dF,F_inverse=F_clamp.Inverse();
    if(!use_constant_ife) dF=F-F_clamp;
    T id_mu_minus_lambda_log_J=id_mu-id_lambda*log(failure_threshold);
    return id_mu*F+id_mu_minus_lambda_log_J*(sqr(F_inverse)*dF-F_inverse)+id_lambda*DIAGONAL_MATRIX<T,d>::Inner_Product(F_inverse,dF)*F_inverse;
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*id_beta*strain_rate+id_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int NEO_HOOKEAN<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void NEO_HOOKEAN<T,d>::
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(id_alpha/TV::dimension+dd*dd)-dd;
    SYMMETRIC_MATRIX<T,d> s=sb*strain_rate+sa*strain_rate.Trace();
    *(MATRIX<T,d>*)aggregate.Get_Array_Pointer()+=s;
}
//#####################################################################
// Function P_From_Strain_Rate_Second_Half
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(id_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,1>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,1>& dP_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    DIAGONAL_MATRIX<T,1> F_inverse=(F.Determinant()>=failure_threshold?F:Clamp_To_Hyperbola(F)).Inverse();
    T mu_minus_lambda_logJ=id_mu+id_lambda*log(F_inverse.Determinant());
    SYMMETRIC_MATRIX<T,1> F_inverse_outer=SYMMETRIC_MATRIX<T,1>::Outer_Product(F_inverse.To_Vector());
    dP_dF.x0000=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x00;//alpha+beta+gamma
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    DIAGONAL_MATRIX<T,2> F_inverse=(F.Determinant()>=failure_threshold?F:Clamp_To_Hyperbola(F)).Inverse();
    T mu_minus_lambda_logJ=id_mu+id_lambda*log(F_inverse.Determinant());
    SYMMETRIC_MATRIX<T,2> F_inverse_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(F_inverse.To_Vector());
    dP_dF.x0000=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x00;//alpha+beta+gamma
    dP_dF.x1111=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;
    dP_dF.x1100=id_lambda*F_inverse_outer.x10;//gamma
    dP_dF.x1010=id_mu;//alpha
    dP_dF.x1001=mu_minus_lambda_logJ*F_inverse_outer.x10;//beta
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    DIAGONAL_MATRIX<T,3> F_inverse=(F.Determinant()>=failure_threshold?F:Clamp_To_Hyperbola(F)).Inverse();
    T mu_minus_lambda_logJ=id_mu+id_lambda*log(F_inverse.Determinant());
    SYMMETRIC_MATRIX<T,3> F_inverse_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(F_inverse.To_Vector());
    dPi_dF.x0000=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x00;
    dPi_dF.x1111=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;
    dPi_dF.x2222=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
    dPi_dF.x1100=id_lambda*F_inverse_outer.x10;
    dPi_dF.x2200=id_lambda*F_inverse_outer.x20;
    dPi_dF.x2211=id_lambda*F_inverse_outer.x21;
    dPi_dF.x1010=id_mu;dPi_dF.x2020=id_mu;dPi_dF.x2121=id_mu;
    dPi_dF.x1001=mu_minus_lambda_logJ*F_inverse_outer.x10;
    dPi_dF.x2002=mu_minus_lambda_logJ*F_inverse_outer.x20;
    dPi_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dPi_dF,const int id) const
{
    Isotropic_Stress_Derivative_Helper(F,dPi_dF,id);
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T NEO_HOOKEAN<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T I1=(F*F.Transposed()).Trace(),J=F.Determinant();
    if(J<=0) return (T)1e16; // TODO: Do something smarter here.
    T log_J=log(J);
    return id_mu*((T).5*(I1-TV::m)-log_J)+(T).5*id_lambda*sqr(log_J);
}
namespace PhysBAM{
template class NEO_HOOKEAN<float,1>;
template class NEO_HOOKEAN<float,2>;
template class NEO_HOOKEAN<float,3>;
template class NEO_HOOKEAN<double,1>;
template class NEO_HOOKEAN<double,2>;
template class NEO_HOOKEAN<double,3>;
}
