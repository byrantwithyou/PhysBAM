//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Alexey Stomakhin, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/cube.h>
#include <Tools/Math_Tools/pow.h>
#include <Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/SVK_EXTRAPOLATED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> SVK_EXTRAPOLATED<T,d>::
SVK_EXTRAPOLATED(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input, const T extra_force_coefficient_input)
    :youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),
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
template<class T,int d> SVK_EXTRAPOLATED<T,d>::
~SVK_EXTRAPOLATED()
{
}
//#####################################################################
// Update Lame Constants
//#####################################################################
template<class T,int d> void SVK_EXTRAPOLATED<T,d>::
Update_Lame_Constants(const T youngs_modulus_input, const T poissons_ratio_input,const T Rayleigh_coefficient_input)
{
    constant_lambda=youngs_modulus_input*poissons_ratio_input/((1+poissons_ratio_input)*(1-2*poissons_ratio_input));
    constant_mu=youngs_modulus_input/(2*(1+poissons_ratio_input));
    constant_alpha=Rayleigh_coefficient_input*constant_lambda;
    constant_beta=Rayleigh_coefficient_input*constant_mu;
    youngs_modulus=youngs_modulus_input; poissons_ratio=poissons_ratio_input;
    
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T SVK_EXTRAPOLATED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    return Energy_Density_Helper(F,simplex);
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T SVK_EXTRAPOLATED<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,2>& F,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return 0;
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T SVK_EXTRAPOLATED<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int simplex) const
{
    T x = F.x11;
    T y = F.x22;
    T z = F.x33;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = constant_mu;
    T la = constant_lambda;   
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;
    
    if ((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        return 1.125*la-0.75*la*sqr(x)-0.75*la*sqr(y)-0.75*la*sqr(z)+0.125*la*(x*x*x*x)+0.25*la*sqr(x)*sqr(y)+0.25*la*sqr(x)*sqr(z)+0.125*la*(y*y*y*y)+0.25*la*sqr(y)*sqr(z)+0.125*la*(z*z*z*z)+0.75*mu-0.5*mu*sqr(x)+0.25*mu*(x*x*x*x)-0.5*mu*sqr(y)+0.25*mu*(y*y*y*y)-0.5*mu*sqr(z)+0.25*mu*(z*z*z*z);
    }
    else if ((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
        return 1.125*la+0.75*mu+0.5*la*a*sqr(y)*x+0.5*la*a*sqr(z)*x+0.125*la*(y*y*y*y)+0.125*la*(z*z*z*z)+0.25*mu*(y*y*y*y)+0.25*mu*(z*z*z*z)+0.75*la*sqr(a)-0.5*mu*sqr(y)-0.75*la*sqr(y)-1.5*la*x*a-0.75*la*sqr(z)-0.5*mu*sqr(z)+0.25*la*sqr(y)*sqr(z)-0.375*la*(a*a*a*a)+0.5*mu*sqr(a)-0.75*mu*(a*a*a*a)+k*sqr(x)+k*sqr(a)+0.5*la*cube(a)*x+mu*cube(a)*x-0.25*la*sqr(a)*sqr(y)-0.25*la*sqr(a)*sqr(z)-2*k*x*a-mu*a*x;
    }
    else if ((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        return -0.5*mu*sqr(x)+1.125*la+0.75*mu+0.5*la*y*a*sqr(x)+0.5*y*la*a*sqr(z)+0.125*la*(x*x*x*x)+0.125*la*(z*z*z*z)+0.25*mu*(x*x*x*x)+0.25*mu*(z*z*z*z)-1.5*la*y*a+0.75*la*sqr(a)-0.75*la*sqr(z)-0.75*la*sqr(x)-0.5*mu*sqr(z)+0.25*la*sqr(x)*sqr(z)-0.375*la*(a*a*a*a)+0.5*mu*sqr(a)-0.75*mu*(a*a*a*a)+k*sqr(a)-0.25*la*sqr(a)*sqr(z)+0.5*y*la*cube(a)-y*mu*a+y*mu*cube(a)-2*y*k*a+k*sqr(y)-0.25*la*sqr(x)*sqr(a);
    }
    else if ((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        return -0.5*mu*sqr(x)+1.125*la+0.75*mu+0.5*la*z*a*sqr(x)+0.5*z*la*a*sqr(y)+0.125*la*(x*x*x*x)+0.125*la*(y*y*y*y)+0.25*mu*(x*x*x*x)+0.25*mu*(y*y*y*y)-1.5*la*z*a+0.75*la*sqr(a)-0.5*mu*sqr(y)-0.75*la*sqr(y)-0.75*la*sqr(x)+0.25*la*sqr(x)*sqr(y)-0.375*la*(a*a*a*a)+0.5*mu*sqr(a)-0.75*mu*(a*a*a*a)+k*sqr(a)-0.25*la*sqr(a)*sqr(y)+0.5*z*la*cube(a)-z*mu*a+z*mu*cube(a)-2*z*k*a+k*sqr(z)-0.25*la*sqr(x)*sqr(a);
    }
    else if ((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        return 1.125*la+0.75*mu+0.5*la*a*sqr(z)*x+x*la*y*sqr(a)+0.5*y*la*a*sqr(z)+0.125*la*(z*z*z*z)+0.25*mu*(z*z*z*z)-1.5*la*y*a+1.5*la*sqr(a)-1.5*la*x*a-0.75*la*sqr(z)-0.5*mu*sqr(z)-0.5*la*(a*a*a*a)+mu*sqr(a)-1.5*mu*(a*a*a*a)+k*sqr(x)+2*k*sqr(a)+mu*cube(a)*x-0.5*la*sqr(a)*sqr(z)-2*k*x*a-mu*a*x-y*mu*a+y*mu*cube(a)-2*y*k*a+k*sqr(y);
    }
    else if ((dx < 0) && (dy >= 0) && (dz < 0)) // Rzx
    {
        return 1.125*la+0.75*mu+0.5*la*a*sqr(y)*x+x*la*z*sqr(a)+0.5*z*la*a*sqr(y)+0.125*la*(y*y*y*y)+0.25*mu*(y*y*y*y)-1.5*la*z*a+1.5*la*sqr(a)-0.5*mu*sqr(y)-0.75*la*sqr(y)-1.5*la*x*a-0.5*la*(a*a*a*a)+mu*sqr(a)-1.5*mu*(a*a*a*a)+k*sqr(x)+2*k*sqr(a)+mu*cube(a)*x-0.5*la*sqr(a)*sqr(y)-2*k*x*a-mu*a*x-z*mu*a+z*mu*cube(a)-2*z*k*a+k*sqr(z);
    }
    else if ((dx >= 0) && (dy < 0) && (dz < 0)) // Ryz
    {
        return -0.5*mu*sqr(x)+1.125*la+0.75*mu+0.5*la*y*a*sqr(x)+0.5*la*z*a*sqr(x)+y*la*z*sqr(a)+0.125*la*(x*x*x*x)+0.25*mu*(x*x*x*x)-1.5*la*z*a-1.5*la*y*a+1.5*la*sqr(a)-0.75*la*sqr(x)-0.5*la*(a*a*a*a)+mu*sqr(a)-1.5*mu*(a*a*a*a)+2*k*sqr(a)-y*mu*a+y*mu*cube(a)-2*y*k*a-z*mu*a+z*mu*cube(a)-2*z*k*a+k*sqr(y)+k*sqr(z)-0.5*la*sqr(x)*sqr(a);
    }
    else // Rxyz
    {
        return 1.125*la+0.75*mu+x*la*y*sqr(a)+x*la*z*sqr(a)+y*la*z*sqr(a)-1.5*la*z*a-1.5*la*y*a+2.25*la*sqr(a)-1.5*la*x*a-0.375*la*(a*a*a*a)+1.5*mu*sqr(a)-2.25*mu*(a*a*a*a)+k*sqr(x)+3*k*sqr(a)-0.5*la*cube(a)*x+mu*cube(a)*x-2*k*x*a-mu*a*x-0.5*y*la*cube(a)-y*mu*a+y*mu*cube(a)-2*y*k*a-0.5*z*la*cube(a)-z*mu*a+z*mu*cube(a)-2*z*k*a+k*sqr(y)+k*sqr(z);
    }
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> SVK_EXTRAPOLATED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    return P_From_Strain_Helper(F,scale,simplex);
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,2> SVK_EXTRAPOLATED<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,2>& F,const T scale,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return DIAGONAL_MATRIX<T,2>();
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> SVK_EXTRAPOLATED<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const
{
    T x = F.x11;
    T y = F.x22;
    T z = F.x33;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = constant_mu;
    T la = constant_lambda;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;
     
    if ((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=-1.5*la*x+0.5*la*cube(x)+0.5*la*x*sqr(y)+0.5*la*x*sqr(z)-mu*x+mu*cube(x);
        result.x22=-1.5*la*y+0.5*la*y*sqr(x)+0.5*la*cube(y)+0.5*la*y*sqr(z)-mu*y+mu*cube(y);
        result.x33=-1.5*la*z+0.5*la*z*sqr(x)+0.5*la*z*sqr(y)+0.5*la*cube(z)-mu*z+mu*cube(z);
        return scale*result;
    }
    else if ((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=-1.5*la*a+0.5*la*cube(a)+0.5*la*a*sqr(y)+0.5*la*a*sqr(z)-mu*a+mu*cube(a)+2*k*x-2*k*a;
        result.x22=-1.5*la*y-0.5*la*y*sqr(a)+0.5*la*cube(y)+0.5*la*y*sqr(z)-mu*y+mu*cube(y)+la*y*a*x;
        result.x33=-1.5*la*z-0.5*la*z*sqr(a)+0.5*la*z*sqr(y)+0.5*la*cube(z)-mu*z+mu*cube(z)+la*z*a*x;
        return scale*result;
    }
    else if ((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=-1.5*la*x+0.5*la*cube(x)-0.5*la*x*sqr(a)+0.5*la*x*sqr(z)-mu*x+mu*cube(x)+la*y*a*x;
        result.x22=-1.5*la*a+0.5*la*a*sqr(x)+0.5*la*cube(a)+0.5*la*a*sqr(z)-mu*a+mu*cube(a)+2*k*y-2*k*a;
        result.x33=-1.5*la*z+0.5*la*z*sqr(x)-0.5*la*z*sqr(a)+0.5*la*cube(z)-mu*z+mu*cube(z)+la*z*a*y;
        return scale*result;
    }
    else if ((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=-1.5*la*x+0.5*la*cube(x)+0.5*la*x*sqr(y)-0.5*la*x*sqr(a)-mu*x+mu*cube(x)+la*z*a*x;
        result.x22=-1.5*la*y+0.5*la*y*sqr(x)+0.5*la*cube(y)-0.5*la*y*sqr(a)-mu*y+mu*cube(y)+la*z*a*y;
        result.x33=-1.5*la*a+0.5*la*a*sqr(x)+0.5*la*a*sqr(y)+0.5*la*cube(a)-mu*a+mu*cube(a)+2*k*z-2*k*a;
        return scale*result;
    }
    else if ((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=-1.5*la*a+0.5*la*a*sqr(z)-mu*a+mu*cube(a)+la*y*sqr(a)+2*k*x-2*k*a;
        result.x22=-1.5*la*a+0.5*la*a*sqr(z)-mu*a+mu*cube(a)+la*x*sqr(a)+2*k*y-2*k*a;
        result.x33=-1.5*la*z-la*z*sqr(a)+0.5*la*cube(z)-mu*z+mu*cube(z)+la*z*a*x+la*z*a*y;
        return scale*result;
    }
    else if ((dx < 0) && (dy >= 0) && (dz < 0)) // Rzx
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=-1.5*la*a+0.5*la*a*sqr(y)-mu*a+mu*cube(a)+la*z*sqr(a)+2*k*x-2*k*a;
        result.x22=-1.5*la*y-la*y*sqr(a)+0.5*la*cube(y)-mu*y+mu*cube(y)+la*z*a*y+la*y*a*x;
        result.x33=-1.5*la*a+0.5*la*a*sqr(y)-mu*a+mu*cube(a)+la*x*sqr(a)+2*k*z-2*k*a;
        return scale*result;
    }
    else if ((dx >= 0) && (dy < 0) && (dz < 0)) // Ryz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=-1.5*la*x+0.5*la*cube(x)-la*x*sqr(a)-mu*x+mu*cube(x)+la*y*a*x+la*z*a*x;
        result.x22=-1.5*la*a+0.5*la*a*sqr(x)-mu*a+mu*cube(a)+la*z*sqr(a)+2*k*y-2*k*a;
        result.x33=-1.5*la*a+0.5*la*a*sqr(x)-mu*a+mu*cube(a)+la*y*sqr(a)+2*k*z-2*k*a;
        return scale*result;
    }
    else // Rxyz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=-1.5*la*a-0.5*la*cube(a)-mu*a+mu*cube(a)+la*y*sqr(a)+la*z*sqr(a)+2*k*x-2*k*a;
        result.x22=-1.5*la*a-0.5*la*cube(a)-mu*a+mu*cube(a)+la*x*sqr(a)+la*z*sqr(a)+2*k*y-2*k*a;
        result.x33=-1.5*la*a-0.5*la*cube(a)-mu*a+mu*cube(a)+la*y*sqr(a)+la*x*sqr(a)+2*k*z-2*k*a;
        return scale*result;
    }
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void SVK_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int triangle) const
{
    return Isotropic_Stress_Derivative_Helper(F,dP_dF,triangle);
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void SVK_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void SVK_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int triangle) const
{
    T x = F.x11;
    T y = F.x22;
    T z = F.x33;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = constant_mu;
    T la = constant_lambda;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;

    T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
    T xmy = x-y; if (fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;
    T xpz = x+z; if (fabs(xpz)<panic_threshold) xpz=xpz<0?-panic_threshold:panic_threshold;
    T xmz = x-z; if (fabs(xmz)<panic_threshold) xmz=xmz<0?-panic_threshold:panic_threshold;
    T ypz = y+z; if (fabs(ypz)<panic_threshold) ypz=ypz<0?-panic_threshold:panic_threshold;
    T ymz = y-z; if (fabs(ymz)<panic_threshold) ymz=ymz<0?-panic_threshold:panic_threshold;

    if ((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        dP_dF.x1111=1.5*la*sqr(x)-1.5*la+0.5*la*sqr(y)+0.5*la*sqr(z)+3*mu*sqr(x)-mu;
        dP_dF.x2211=la*y*x;
        dP_dF.x2222=1.5*la*sqr(y)-1.5*la+0.5*la*sqr(x)+0.5*la*sqr(z)+3*mu*sqr(y)-mu;
        dP_dF.x3311=la*x*z;
        dP_dF.x3322=la*z*y;
        dP_dF.x3333=1.5*la*sqr(z)-1.5*la+0.5*la*sqr(x)+0.5*la*sqr(y)+3*mu*sqr(z)-mu;
        dP_dF.x2112=x*mu*y;
        dP_dF.x2121=0.5*la*sqr(x)+mu*sqr(x)-1.5*la+0.5*la*sqr(z)-mu+0.5*la*sqr(y)+mu*sqr(y);
        dP_dF.x3113=x*mu*z;
        dP_dF.x3131=0.5*la*sqr(x)+mu*sqr(x)-1.5*la+0.5*la*sqr(y)-mu+0.5*la*sqr(z)+mu*sqr(z);
        dP_dF.x3223=y*mu*z;
        dP_dF.x3232=0.5*la*sqr(y)+mu*sqr(y)-1.5*la+0.5*la*sqr(x)-mu+0.5*la*sqr(z)+mu*sqr(z);
    }
    else if ((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
        dP_dF.x1111=2*k;
        dP_dF.x2211=la*y*a;
        dP_dF.x2222=1.5*la*sqr(y)-1.5*la-0.5*la*sqr(a)+0.5*la*sqr(z)+3*mu*sqr(y)-mu+la*x*a;
        dP_dF.x3311=la*z*a;
        dP_dF.x3322=la*z*y;
        dP_dF.x3333=1.5*la*sqr(z)-1.5*la-0.5*la*sqr(a)+0.5*la*sqr(y)+3*mu*sqr(z)-mu+la*x*a;
        dP_dF.x2112=0.5*y*(3*la*x+la*x*sqr(a)-la*x*sqr(y)-la*x*sqr(z)+2*mu*x-2*x*mu*sqr(y)-2*la*a*sqr(x)-3*la*a+la*cube(a)+la*a*sqr(y)+la*a*sqr(z)-2*mu*a+2*mu*cube(a)+4*k*x-4*k*a)/(xpy*xmy);
        dP_dF.x2121=0.5*(3*la*sqr(y)+la*sqr(a)*sqr(y)-la*(y*y*y*y)-la*sqr(y)*sqr(z)+2*mu*sqr(y)-2*mu*(y*y*y*y)-la*a*sqr(y)*x-3*la*x*a+la*cube(a)*x+la*a*sqr(z)*x-2*mu*a*x+2*mu*cube(a)*x+4*k*sqr(x)-4*k*x*a)/(xpy*xmy);
        dP_dF.x3113=0.5*z*(3*la*x+la*x*sqr(a)-la*x*sqr(y)-la*x*sqr(z)+2*mu*x-2*x*mu*sqr(z)-2*la*a*sqr(x)-3*la*a+la*cube(a)+la*a*sqr(y)+la*a*sqr(z)-2*mu*a+2*mu*cube(a)+4*k*x-4*k*a)/(xpz*xmz);
        dP_dF.x3131=0.5*(3*la*sqr(z)+la*sqr(a)*sqr(z)-la*sqr(y)*sqr(z)-la*(z*z*z*z)+2*mu*sqr(z)-2*mu*(z*z*z*z)-la*a*sqr(z)*x-3*la*x*a+la*cube(a)*x+la*a*sqr(y)*x-2*mu*a*x+2*mu*cube(a)*x+4*k*sqr(x)-4*k*x*a)/(xpz*xmz);
        dP_dF.x3223=y*mu*z;
        dP_dF.x3232=0.5*la*sqr(y)+mu*sqr(y)-1.5*la-0.5*la*sqr(a)-mu+la*x*a+0.5*la*sqr(z)+mu*sqr(z);
    }
    else if ((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        dP_dF.x1111=1.5*la*sqr(x)-1.5*la-0.5*la*sqr(a)+0.5*la*sqr(z)+3*mu*sqr(x)-mu+la*y*a;
        dP_dF.x2211=la*x*a;
        dP_dF.x2222=2*k;
        dP_dF.x3311=la*x*z;
        dP_dF.x3322=la*z*a;
        dP_dF.x3333=1.5*la*sqr(z)-1.5*la+0.5*la*sqr(x)-0.5*la*sqr(a)+3*mu*sqr(z)-mu+la*y*a;
        dP_dF.x2112=-0.5*x*(-3*la*a+la*a*sqr(x)+la*cube(a)+la*a*sqr(z)-2*mu*a+2*mu*cube(a)+4*k*y-4*k*a+3*la*y-la*y*sqr(x)+la*y*sqr(a)-la*y*sqr(z)+2*mu*y-2*y*mu*sqr(x)-2*la*a*sqr(y))/(xpy*xmy);
        dP_dF.x2121=-0.5*(-3*la*y*a-la*y*a*sqr(x)+y*la*cube(a)+y*la*a*sqr(z)-2*y*mu*a+2*y*mu*cube(a)+4*k*sqr(y)-4*y*k*a+3*la*sqr(x)-la*(x*x*x*x)+la*sqr(x)*sqr(a)-la*sqr(x)*sqr(z)+2*mu*sqr(x)-2*mu*(x*x*x*x))/(xpy*xmy);
        dP_dF.x3113=x*mu*z;
        dP_dF.x3131=0.5*la*sqr(x)+mu*sqr(x)-1.5*la+la*y*a-0.5*la*sqr(a)-mu+0.5*la*sqr(z)+mu*sqr(z);
        dP_dF.x3223=0.5*z*(3*la*y-la*y*sqr(x)+la*y*sqr(a)-la*y*sqr(z)+2*mu*y-2*y*mu*sqr(z)-2*la*a*sqr(y)-3*la*a+la*a*sqr(x)+la*cube(a)+la*a*sqr(z)-2*mu*a+2*mu*cube(a)+4*k*y-4*k*a)/(ypz*ymz);
        dP_dF.x3232=0.5*(3*la*sqr(z)-la*sqr(x)*sqr(z)+la*sqr(a)*sqr(z)-la*(z*z*z*z)+2*mu*sqr(z)-2*mu*(z*z*z*z)-y*la*a*sqr(z)-3*la*y*a+la*y*a*sqr(x)+y*la*cube(a)-2*y*mu*a+2*y*mu*cube(a)+4*k*sqr(y)-4*y*k*a)/(ypz*ymz);
    }
    else if ((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        dP_dF.x1111=1.5*la*sqr(x)-1.5*la+0.5*la*sqr(y)-0.5*la*sqr(a)+3*mu*sqr(x)-mu+la*z*a;
        dP_dF.x2211=la*y*x;
        dP_dF.x2222=1.5*la*sqr(y)-1.5*la+0.5*la*sqr(x)-0.5*la*sqr(a)+3*mu*sqr(y)-mu+la*z*a;
        dP_dF.x3311=la*x*a;
        dP_dF.x3322=la*y*a;
        dP_dF.x3333=2*k;
        dP_dF.x2112=x*mu*y;
        dP_dF.x2121=0.5*la*sqr(x)+mu*sqr(x)-1.5*la+la*z*a-0.5*la*sqr(a)-mu+0.5*la*sqr(y)+mu*sqr(y);
        dP_dF.x3113=-0.5*x*(-3*la*a+la*a*sqr(x)+la*a*sqr(y)+la*cube(a)-2*mu*a+2*mu*cube(a)+4*k*z-4*k*a+3*la*z-la*z*sqr(x)-la*z*sqr(y)+la*z*sqr(a)+2*mu*z-2*z*mu*sqr(x)-2*la*a*sqr(z))/(xpz*xmz);
        dP_dF.x3131=-0.5*(-3*la*z*a-la*z*a*sqr(x)+z*la*a*sqr(y)+z*la*cube(a)-2*z*mu*a+2*z*mu*cube(a)+4*k*sqr(z)-4*z*k*a+3*la*sqr(x)-la*(x*x*x*x)-la*sqr(x)*sqr(y)+la*sqr(x)*sqr(a)+2*mu*sqr(x)-2*mu*(x*x*x*x))/(xpz*xmz);
        dP_dF.x3223=-0.5*y*(-3*la*a+la*a*sqr(x)+la*a*sqr(y)+la*cube(a)-2*mu*a+2*mu*cube(a)+4*k*z-4*k*a+3*la*z-la*z*sqr(x)-la*z*sqr(y)+la*z*sqr(a)+2*mu*z-2*z*mu*sqr(y)-2*la*a*sqr(z))/(ypz*ymz);
        dP_dF.x3232=-0.5*(-3*la*z*a+la*z*a*sqr(x)-z*la*a*sqr(y)+z*la*cube(a)-2*z*mu*a+2*z*mu*cube(a)+4*k*sqr(z)-4*z*k*a+3*la*sqr(y)-la*sqr(x)*sqr(y)-la*(y*y*y*y)+la*sqr(a)*sqr(y)+2*mu*sqr(y)-2*mu*(y*y*y*y))/(ypz*ymz);
    }
    else if ((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        dP_dF.x1111=2*k;
        dP_dF.x2211=la*sqr(a);
        dP_dF.x2222=2*k;
        dP_dF.x3311=la*z*a;
        dP_dF.x3322=la*z*a;
        dP_dF.x3333=1.5*la*sqr(z)-1.5*la-la*sqr(a)+3*mu*sqr(z)-mu+la*x*a+la*y*a;
        dP_dF.x2112=-0.5*(2*la*x*a-3*la+la*sqr(z)-4*k-2*mu+2*mu*sqr(a)+2*la*y*a)*a/xpy;
        dP_dF.x2121=0.5*(4*k*x-3*la*a+la*a*sqr(z)-2*mu*a+2*mu*cube(a)-4*k*a+4*k*y)/xpy;
        dP_dF.x3113=0.5*z*(3*la*x+2*la*x*sqr(a)-la*x*sqr(z)+2*mu*x-2*x*mu*sqr(z)-2*la*a*sqr(x)-2*la*y*a*x-3*la*a+la*a*sqr(z)-2*mu*a+2*mu*cube(a)+2*la*y*sqr(a)+4*k*x-4*k*a)/(xpz*xmz);
        dP_dF.x3131=0.5*(3*la*sqr(z)+2*la*sqr(a)*sqr(z)-la*(z*z*z*z)+2*mu*sqr(z)-2*mu*(z*z*z*z)-la*a*sqr(z)*x-2*y*la*a*sqr(z)-3*la*x*a-2*mu*a*x+2*mu*cube(a)*x+2*x*la*y*sqr(a)+4*k*sqr(x)-4*k*x*a)/(xpz*xmz);
        dP_dF.x3223=0.5*z*(3*la*y+2*la*y*sqr(a)-la*y*sqr(z)+2*mu*y-2*y*mu*sqr(z)-2*la*y*a*x-2*la*a*sqr(y)-3*la*a+la*a*sqr(z)-2*mu*a+2*mu*cube(a)+2*la*x*sqr(a)+4*k*y-4*k*a)/(ypz*ymz);
        dP_dF.x3232=0.5*(3*la*sqr(z)+2*la*sqr(a)*sqr(z)-la*(z*z*z*z)+2*mu*sqr(z)-2*mu*(z*z*z*z)-2*la*a*sqr(z)*x-y*la*a*sqr(z)-3*la*y*a-2*y*mu*a+2*y*mu*cube(a)+2*x*la*y*sqr(a)+4*k*sqr(y)-4*y*k*a)/(ypz*ymz);
    }
    else if ((dx < 0) && (dy >= 0) && (dz < 0)) // Rzx
    {
        dP_dF.x1111=2*k;
        dP_dF.x2211=la*y*a;
        dP_dF.x2222=1.5*la*sqr(y)-1.5*la-la*sqr(a)+3*mu*sqr(y)-mu+la*z*a+la*x*a;
        dP_dF.x3311=la*sqr(a);
        dP_dF.x3322=la*y*a;
        dP_dF.x3333=2*k;
        dP_dF.x2112=0.5*y*(3*la*x+2*la*x*sqr(a)-la*x*sqr(y)+2*mu*x-2*x*mu*sqr(y)-2*la*z*a*x-2*la*a*sqr(x)-3*la*a+la*a*sqr(y)-2*mu*a+2*mu*cube(a)+2*la*z*sqr(a)+4*k*x-4*k*a)/(xpy*xmy);
        dP_dF.x2121=0.5*(3*la*sqr(y)+2*la*sqr(a)*sqr(y)-la*(y*y*y*y)+2*mu*sqr(y)-2*mu*(y*y*y*y)-2*z*la*a*sqr(y)-la*a*sqr(y)*x-3*la*x*a-2*mu*a*x+2*mu*cube(a)*x+2*x*la*z*sqr(a)+4*k*sqr(x)-4*k*x*a)/(xpy*xmy);
        dP_dF.x3113=-0.5*(2*la*x*a-3*la+la*sqr(y)-4*k-2*mu+2*mu*sqr(a)+2*la*z*a)*a/xpz;
        dP_dF.x3131=0.5*(4*k*x-3*la*a+la*a*sqr(y)-2*mu*a+2*mu*cube(a)-4*k*a+4*k*z)/xpz;
        dP_dF.x3223=-0.5*y*(-3*la*a+la*a*sqr(y)-2*mu*a+2*mu*cube(a)+2*la*x*sqr(a)+4*k*z-4*k*a+3*la*z+2*la*z*sqr(a)-la*z*sqr(y)+2*mu*z-2*z*mu*sqr(y)-2*la*a*sqr(z)-2*la*z*a*x)/(ypz*ymz);
        dP_dF.x3232=-0.5*(-3*la*z*a-z*la*a*sqr(y)-2*z*mu*a+2*z*mu*cube(a)+2*x*la*z*sqr(a)+4*k*sqr(z)-4*z*k*a+3*la*sqr(y)+2*la*sqr(a)*sqr(y)-la*(y*y*y*y)+2*mu*sqr(y)-2*mu*(y*y*y*y)-2*la*a*sqr(y)*x)/(ypz*ymz);
    }
    else if ((dx >= 0) && (dy < 0) && (dz < 0)) // Ryz
    {
        dP_dF.x1111=1.5*la*sqr(x)-1.5*la-la*sqr(a)+3*mu*sqr(x)-mu+la*y*a+la*z*a;
        dP_dF.x2211=la*x*a;
        dP_dF.x2222=2*k;
        dP_dF.x3311=la*x*a;
        dP_dF.x3322=la*sqr(a);
        dP_dF.x3333=2*k;
        dP_dF.x2112=-0.5*x*(-3*la*a+la*a*sqr(x)-2*mu*a+2*mu*cube(a)+2*la*z*sqr(a)+4*k*y-4*k*a+3*la*y-la*y*sqr(x)+2*la*y*sqr(a)+2*mu*y-2*y*mu*sqr(x)-2*la*a*sqr(y)-2*la*z*a*y)/(xpy*xmy);
        dP_dF.x2121=-0.5*(-3*la*y*a-la*y*a*sqr(x)-2*y*mu*a+2*y*mu*cube(a)+2*y*la*z*sqr(a)+4*k*sqr(y)-4*y*k*a+3*la*sqr(x)-la*(x*x*x*x)+2*la*sqr(x)*sqr(a)+2*mu*sqr(x)-2*mu*(x*x*x*x)-2*la*z*a*sqr(x))/(xpy*xmy);
        dP_dF.x3113=-0.5*x*(-3*la*a+la*a*sqr(x)-2*mu*a+2*mu*cube(a)+2*la*y*sqr(a)+4*k*z-4*k*a+3*la*z-la*z*sqr(x)+2*la*z*sqr(a)+2*mu*z-2*z*mu*sqr(x)-2*la*z*a*y-2*la*a*sqr(z))/(xpz*xmz);
        dP_dF.x3131=-0.5*(-3*la*z*a-la*z*a*sqr(x)-2*z*mu*a+2*z*mu*cube(a)+2*y*la*z*sqr(a)+4*k*sqr(z)-4*z*k*a+3*la*sqr(x)-la*(x*x*x*x)+2*la*sqr(x)*sqr(a)+2*mu*sqr(x)-2*mu*(x*x*x*x)-2*la*y*a*sqr(x))/(xpz*xmz);
        dP_dF.x3223=-0.5*(2*la*y*a-3*la+la*sqr(x)-4*k-2*mu+2*mu*sqr(a)+2*la*z*a)*a/ypz;
        dP_dF.x3232=0.5*(4*k*y-3*la*a+la*a*sqr(x)-2*mu*a+2*mu*cube(a)-4*k*a+4*k*z)/ypz;
    }
    else // Rxyz
    {
        dP_dF.x1111=2*k;
        dP_dF.x2211=la*sqr(a);
        dP_dF.x2222=2*k;
        dP_dF.x3311=la*sqr(a);
        dP_dF.x3322=la*sqr(a);
        dP_dF.x3333=2*k;
        dP_dF.x2112=0.5*(-2*la*x*a+la*sqr(a)+2*mu-2*la*z*a+4*k-2*mu*sqr(a)+3*la-2*la*y*a)*a/xpy;
        dP_dF.x2121=-0.5*(-4*k*x+3*la*a+la*cube(a)+2*mu*a-2*mu*cube(a)-2*la*z*sqr(a)+4*k*a-4*k*y)/xpy;
        dP_dF.x3113=0.5*(-2*la*x*a+la*sqr(a)+2*mu-2*la*z*a+4*k-2*mu*sqr(a)+3*la-2*la*y*a)*a/xpz;
        dP_dF.x3131=-0.5*(-4*k*x+3*la*a+la*cube(a)+2*mu*a-2*mu*cube(a)-2*la*y*sqr(a)+4*k*a-4*k*z)/xpz;
        dP_dF.x3223=0.5*(-2*la*x*a+la*sqr(a)+2*mu-2*la*z*a+4*k-2*mu*sqr(a)+3*la-2*la*y*a)*a/ypz;
        dP_dF.x3232=-0.5*(-4*k*y+3*la*a+la*cube(a)+2*mu*a-2*mu*cube(a)-2*la*x*sqr(a)+4*k*a-4*k*z)/ypz;
    }
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> SVK_EXTRAPOLATED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int SVK_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void SVK_EXTRAPOLATED<T,d>::
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
template<class T,int d> MATRIX<T,d> SVK_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
namespace PhysBAM{
template class SVK_EXTRAPOLATED<float,2>;
template class SVK_EXTRAPOLATED<float,3>;
template class SVK_EXTRAPOLATED<double,2>;
template class SVK_EXTRAPOLATED<double,3>;
}
