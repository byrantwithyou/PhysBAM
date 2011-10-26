//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Alexey Stomakhin, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_COROTATED_BLEND.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
NEO_HOOKEAN_COROTATED_BLEND(const T neo_youngs_modulus_input,
                         const T neo_poissons_ratio_input,
                         const T neo_Rayleigh_coefficient):
    neo_youngs_modulus(neo_youngs_modulus_input),neo_poissons_ratio(neo_poissons_ratio_input),
    panic_threshold((T)1e-6)
{
    assert(neo_poissons_ratio>-1&&neo_poissons_ratio<.5);
    
    constant_lambda = neo_youngs_modulus*neo_poissons_ratio/((1+neo_poissons_ratio)*(1-2*neo_poissons_ratio));
    constant_mu     = neo_youngs_modulus/(2*(1+neo_poissons_ratio));
    constant_alpha  = neo_Rayleigh_coefficient*constant_lambda;
    constant_beta   = neo_Rayleigh_coefficient*constant_mu;

    neo_mu     = constant_mu;
    neo_lambda = constant_lambda;

    /*** note ***/ // TODO: Should choose the constants here based on the input parameters 

        cor_mu     = constant_mu*2;
        cor_lambda = constant_lambda*2;

        J_min = 0.3;
        J_max = 0.9;

    /*** end note ***/

    heaviside.Initialize(J_min,J_max);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
~NEO_HOOKEAN_COROTATED_BLEND()
{
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    return Energy_Density_Helper(F,simplex);
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,2>& F,const int simplex) const
{
    T J=F.Determinant();

    if (J>=J_max) // Neo Hookean
    {
        T I1=(F*F.Transposed()).Trace();
        T log_J=log(J);
        return neo_mu*((T).5*(I1-TV::m)-log_J)+(T).5*neo_lambda*sqr(log_J);
    }
    else if (J<=J_min) // Corotated
    {
        DIAGONAL_MATRIX<T,2> Fm1=F-1;
        return cor_mu*(Fm1*Fm1).Trace()+(T).5*cor_lambda*sqr(Fm1.Trace());
    }
    else // Transition Region J_min < J < J_max
    {
        T I1=(F*F.Transposed()).Trace();
        T log_J=log(J);
        T neo_energy = neo_mu*((T).5*(I1-TV::m)-log_J)+(T).5*neo_lambda*sqr(log_J);
        
        DIAGONAL_MATRIX<T,2> Fm1=F-1;
        T cor_energy = cor_mu*(Fm1*Fm1).Trace()+(T).5*cor_lambda*sqr(Fm1.Trace());
        
        T t = heaviside.H(J);

        return t*neo_energy + (1-t)*cor_energy;
    }
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return 0;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    return P_From_Strain_Helper(F,scale,simplex);
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,2> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,2>& F,const T scale,const int simplex) const
{
    T J=F.Determinant();

    if (J>=J_max) // Neo Hookean
    {  
        return scale*(neo_mu*F-(neo_mu-neo_lambda*log(J))*F.Inverse());
    }
    else if (J<=J_min) // Corotated
    {
        DIAGONAL_MATRIX<T,2> Fm1=F-1;
        return scale*(2*cor_mu*Fm1+cor_lambda*Fm1.Trace());
    }
    else // Transition Region J_min < J < J_max
    {
        T I1=(F*F.Transposed()).Trace();
        T log_J=log(J);
        T neo_energy = neo_mu*((T).5*(I1-TV::m)-log_J)+(T).5*neo_lambda*sqr(log_J);
        DIAGONAL_MATRIX<T,2> grad_neo_energy = neo_mu*F-(neo_mu-neo_lambda*log_J)*F.Inverse();

        DIAGONAL_MATRIX<T,2> Fm1=F-1;
        T cor_energy = cor_mu*(Fm1*Fm1).Trace()+(T).5*cor_lambda*sqr(Fm1.Trace());
        DIAGONAL_MATRIX<T,2> grad_cor_energy = 2*cor_mu*Fm1+cor_lambda*Fm1.Trace();
        
        T t = heaviside.H(J);
        T tJ = heaviside.HJ(J);
        DIAGONAL_MATRIX<T,2> grad_t;
        grad_t.x11 = tJ*F.x22;
        grad_t.x22 = tJ*F.x11;

        return scale*(grad_neo_energy*t + grad_cor_energy*(1-t) + grad_t*(neo_energy-cor_energy));
    }
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return DIAGONAL_MATRIX<T,3>();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int triangle) const
{
    return Isotropic_Stress_Derivative_Helper(F,dP_dF,triangle);
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    T J=F.Determinant();

    if (J>=J_max) // Neo Hookean
    {     
        DIAGONAL_MATRIX<T,2> F_inverse=F.Inverse();
        T mu_minus_lambda_logJ=neo_mu+neo_lambda*log(F_inverse.Determinant());
        SYMMETRIC_MATRIX<T,2> F_inverse_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(F_inverse.To_Vector());
        dP_dF.x1111=neo_mu+(neo_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;//alpha+beta+gamma
        dP_dF.x2222=neo_mu+(neo_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
        dP_dF.x2211=neo_lambda*F_inverse_outer.x21;//gamma
        dP_dF.x2121=neo_mu;//alpha
        dP_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;//beta
    }
    else if (J<=J_min) // Corotated
    {
         T mu=cor_mu,la=cor_lambda,mu2la=2*mu+la,la2mu2=2*la+2*mu;
         T d12=F.x11+F.x22;if(fabs(d12)<panic_threshold) d12=d12<0?-panic_threshold:panic_threshold;
         T i12=la2mu2/d12;
         dP_dF.x1111=mu2la;
         dP_dF.x2112=i12-la;
         dP_dF.x2121=mu2la-i12;
         dP_dF.x2211=la;
         dP_dF.x2222=mu2la;
    }
    else // Transition Region J_min < J < J_max
    {

    }

    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int triangle) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    SYMMETRIC_MATRIX<T,d> s=sb*strain_rate+sa*strain_rate.Trace();
    *(MATRIX<T,d>*)aggregate.Get_Array_Pointer()+=s;
}
//#####################################################################
// Function P_From_Strain_Rate_Second_Half
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
template class NEO_HOOKEAN_COROTATED_BLEND<float,2>;
template class NEO_HOOKEAN_COROTATED_BLEND<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NEO_HOOKEAN_COROTATED_BLEND<double,2>;
template class NEO_HOOKEAN_COROTATED_BLEND<double,3>;
#endif
