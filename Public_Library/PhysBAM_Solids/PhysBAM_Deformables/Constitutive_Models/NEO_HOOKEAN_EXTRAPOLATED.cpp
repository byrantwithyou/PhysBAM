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
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED<T,d>::
NEO_HOOKEAN_EXTRAPOLATED(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input, const T extra_force_coefficient_input)
    :youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),extrapolation_cutoff(extrapolation_cutoff_input),extra_force_coefficient(extra_force_coefficient_input)
{
    assert(poissons_ratio>-1&&poissons_ratio<.5);
    constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    constant_mu=youngs_modulus/(2*(1+poissons_ratio));
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;
    base.Initialize(constant_mu,constant_lambda);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED<T,d>::
~NEO_HOOKEAN_EXTRAPOLATED()
{
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    if (d==2)
    {
        T x = F.x1;
        T y = F.x2;

        T dx = x - extrapolation_cutoff;
        T dy = y - extrapolation_cutoff;

        if ((dx >= 0) && (dy >= 0))
        {
            T I1=(F*F.Transposed()).Trace(),J=F.Determinant();
            T log_J=log(J);
            return constant_mu*((T).5*(I1-TV::m)-log_J)+(T).5*constant_lambda*sqr(log_J);
        }
        else if ((dx < 0) && (dy >= 0))
        {
            return 
        }
        else if ((dx >= 0) && (dy < 0))
        {

        }
        else // ((dx < 0) && (dy < 0))
        {

        }
    }
    else PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
// clamp to hyperbola to avoid indefiniteness "automatically"
template<class T,int d> DIAGONAL_MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda,J=F.Determinant();
    if (J>=failure_threshold) return scale_mu*F-(scale_mu-scale_lambda*log(J))*F.Inverse();
    return F;
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int NEO_HOOKEAN_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED<T,d>::
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
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    DIAGONAL_MATRIX<T,2> F_inverse=F.Clamp_Min(failure_threshold).Inverse();
    T mu_minus_lambda_logJ=constant_mu+constant_lambda*log(F_inverse.Determinant());
    SYMMETRIC_MATRIX<T,2> F_inverse_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(F_inverse.To_Vector());
    dP_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;//alpha+beta+gamma
    dP_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
    dP_dF.x2211=constant_lambda*F_inverse_outer.x21;//gamma
    dP_dF.x2121=constant_mu;//alpha
    dP_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;//beta
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const
{
    DIAGONAL_MATRIX<T,3> F_inverse=F.Clamp_Min(failure_threshold).Inverse();
    T mu_minus_lambda_logJ=constant_mu+constant_lambda*log(F_inverse.Determinant());
    SYMMETRIC_MATRIX<T,3> F_inverse_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(F_inverse.To_Vector());
    dPi_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;
    dPi_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
    dPi_dF.x3333=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x33;
    dPi_dF.x2211=constant_lambda*F_inverse_outer.x21;
    dPi_dF.x3311=constant_lambda*F_inverse_outer.x31;
    dPi_dF.x3322=constant_lambda*F_inverse_outer.x32;
    dPi_dF.x2121=constant_mu;dPi_dF.x3131=constant_mu;dPi_dF.x3232=constant_mu;
    dPi_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;
    dPi_dF.x3113=mu_minus_lambda_logJ*F_inverse_outer.x31;
    dPi_dF.x3223=mu_minus_lambda_logJ*F_inverse_outer.x32;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
}
template class NEO_HOOKEAN_EXTRAPOLATED<float,2>;
template class NEO_HOOKEAN_EXTRAPOLATED<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NEO_HOOKEAN_EXTRAPOLATED<double,2>;
template class NEO_HOOKEAN_EXTRAPOLATED<double,3>;
#endif
