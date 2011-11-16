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
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_SMOOTH.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
NEO_HOOKEAN_EXTRAPOLATED_SMOOTH(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input, const T extra_force_coefficient_input):
    youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),
    extrapolation_cutoff(extrapolation_cutoff_input),extra_force_coefficient(extra_force_coefficient_input),
    panic_threshold((T)1e-6)
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
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
~NEO_HOOKEAN_EXTRAPOLATED_SMOOTH()
{
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    return Energy_Density_Helper(F,simplex);
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,2>& F,const int simplex) const
{
    T x = F.x11;
    T y = F.x22;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient;
    
    if ((dx >= 0) && (dy >= 0))
    {
        T I1=(F*F.Transposed()).Trace(),J=F.Determinant();
        T log_J=log(J);
        return constant_mu*((T).5*(I1-TV::m)-log_J)+(T).5*constant_lambda*sqr(log_J);
    }
    else if ((dx < 0) && (dy >= 0))
    {
        return base.E(a,y) + base.Ex(a,y)*dx + k*dx*dx;
    }
    else if ((dx >= 0) && (dy < 0))
    {
        return base.E(x,a) + base.Ey(x,a)*dy + k*dy*dy;
    }
    else // ((dx < 0) && (dy < 0))
    {
        return base.E(a,a) + base.Ex(a,a)*dx + base.Ey(a,a)*dy + base.Exy(a,a)*dx*dy + k*dx*dx + k*dy*dy;
    }
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return 0;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    return P_From_Strain_Helper(F,scale,simplex);
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,2> NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,2>& F,const T scale,const int simplex) const
{
    T x = F.x11;
    T y = F.x22;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient;
    
    if ((dx >= 0) && (dy >= 0))
    {  
        T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda,J=F.Determinant();
        return scale_mu*F-(scale_mu-scale_lambda*log(J))*F.Inverse();
    }
    else if ((dx < 0) && (dy >= 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = base.Ex(a,y) + 2*k*dx;
        result.x22 = base.Ey(a,y) + base.Exy(a,y)*dx;
        return scale*result;
    }
    else if ((dx >= 0) && (dy < 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = base.Ex(x,a) + base.Exy(x,a)*dy;
        result.x22 = base.Ey(x,a) + 2*k*dy;
        return scale*result;
    }
    else // ((dx < 0) && (dy < 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = base.Ex(a,a) + base.Exy(a,a)*dy + 2*k*dx;
        result.x22 = base.Ey(a,a) + base.Exy(a,a)*dx + 2*k*dy;
        return scale*result;
    }
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return DIAGONAL_MATRIX<T,3>();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int triangle) const
{
    return Isotropic_Stress_Derivative_Helper(F,dP_dF,triangle);
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    T x = F.x11;
    T y = F.x22;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient;
    
    if ((dx >= 0) && (dy >= 0))
    {     
        DIAGONAL_MATRIX<T,2> F_inverse=F.Inverse();
        T mu_minus_lambda_logJ=constant_mu+constant_lambda*log(F_inverse.Determinant());
        SYMMETRIC_MATRIX<T,2> F_inverse_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(F_inverse.To_Vector());
        dP_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;//alpha+beta+gamma
        dP_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
        dP_dF.x2211=constant_lambda*F_inverse_outer.x21;//gamma
        dP_dF.x2121=constant_mu;//alpha
        dP_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;//beta
    }
    else if ((dx < 0) && (dy >= 0))
    {
        dP_dF.x1111=2*k;
        dP_dF.x2222=base.Eyy(a,y)+base.Exyy(a,y)*dx;
        dP_dF.x2211=base.Exy(a,y);
        
        T Ex = base.Ex(a,y)+2*k*dx;
        T Ey = base.Ey(a,y)+base.Exy(a,y)*dx;

        T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T xmy = x-y; if (fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;

        dP_dF.x2112=(-Ey*x+Ex*y)/(xpy*xmy);
        dP_dF.x2121=(-Ey*y+Ex*x)/(xpy*xmy);
    }
    else if ((dx >= 0) && (dy < 0))
    {
        dP_dF.x1111=base.Exx(x,a)+base.Exxy(x,a)*dy;
        dP_dF.x2222=2*k;
        dP_dF.x2211=base.Exy(x,a);
        
        T Ex = base.Ex(x,a)+base.Exy(x,a)*dy;
        T Ey = base.Ey(x,a)+2*k*dy;

        T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T xmy = x-y; if (fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;

        dP_dF.x2112=(-Ey*x+Ex*y)/(xpy*xmy);
        dP_dF.x2121=(-Ey*y+Ex*x)/(xpy*xmy);
    }
    else // ((dx < 0) && (dy < 0))
    {
        dP_dF.x1111=2*k;
        dP_dF.x2222=2*k;
        dP_dF.x2211=base.Exy(a,a);

        T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T cmn = (base.Ex(a,a)-base.Exy(a,a)*a-2*k*a)/xpy;

        dP_dF.x2112=-cmn-base.Exy(a,a);
        dP_dF.x2121=cmn+2*k;
    }

    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int triangle) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
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
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
template class NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<float,2>;
template class NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<double,2>;
template class NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<double,3>;
#endif
