//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Alexey Stomakhin, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
// #####################################################################
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
NEO_HOOKEAN_EXTRAPOLATED_SMOOTH(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input):
    youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),
    extrapolation_cutoff(extrapolation_cutoff_input),
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

    T mu = constant_mu;
    T la = constant_lambda;
    
    T a = extrapolation_cutoff;
     
    if ((dx >= 0) && (dy >= 0))
    {
        T I1=(F*F.Transposed()).Trace(),J=F.Determinant();
        T log_J=log(J);
        return constant_mu*((T).5*(I1-TV::m)-log_J)+(T).5*constant_lambda*sqr(log_J);
    }
    else if ((dx < 0) && (dy >= 0))
    {
        return 0.5*(mu*sqr(y)*sqr(a)+mu*sqr(a)-2*mu*sqr(a)*log(y*a)+la*sqr(a)*sqr(log(y*a))-3*sqr(a)*la*log(y*a)-4*mu*x*a+4*a*la*log(y*a)*x+la*sqr(a)-2*x*a*la+sqr(x)*mu*sqr(a)+mu*sqr(x)+sqr(x)*la-la*log(y*a)*sqr(x))/sqr(a);
    }
    else if ((dx >= 0) && (dy < 0))
    {
        return 0.5*(sqr(x)*mu*sqr(a)+mu*sqr(a)-2*mu*sqr(a)*log(x*a)+la*sqr(a)*sqr(log(x*a))-3*sqr(a)*la*log(x*a)-4*mu*y*a+4*a*y*la*log(x*a)+la*sqr(a)-2*y*a*la+mu*sqr(y)*sqr(a)+mu*sqr(y)+sqr(y)*la-sqr(y)*la*log(x*a))/sqr(a);
    }
    else // ((dx < 0) && (dy < 0))
    {
        return 0.25*(2*mu*sqr(y)*sqr(a)-8*mu*y*cube(a)-8*x*mu*cube(a)+2*sqr(x)*mu*sqr(a)+8*mu*(a*a*a*a)+8*cube(a)*y*la*log(sqr(a))-2*sqr(a)*sqr(y)*la*log(sqr(a))+8*x*cube(a)*la*log(sqr(a))+16*x*y*sqr(a)*la-4*a*sqr(y)*la*x-2*sqr(x)*sqr(a)*la*log(sqr(a))-4*sqr(x)*y*a*la+13*(a*a*a*a)*la+5*sqr(a)*sqr(y)*la+2*(a*a*a*a)*mu*sqr(y)-16*cube(a)*y*la-12*(a*a*a*a)*la*log(sqr(a))-4*(a*a*a*a)*mu*log(sqr(a))+2*(a*a*a*a)*la*sqr(log(sqr(a)))-16*x*la*cube(a)+5*sqr(x)*la*sqr(a)+2*sqr(x)*mu*(a*a*a*a)+sqr(x)*sqr(y)*la)/(a*a*a*a);
    }
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int simplex) const
{
    T x = F.x11;
    T y = F.x22;
    T z = F.x22;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = constant_mu;
    T la = constant_lambda;
    
    T a = extrapolation_cutoff;
     
    if ((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        T I1=(F*F.Transposed()).Trace(),J=F.Determinant();
        T log_J=log(J);
        return constant_mu*((T).5*(I1-TV::m)-log_J)+(T).5*constant_lambda*sqr(log_J);
    }
    else if ((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
        return 1/2*(mu*y^2*a^2+mu*z^2*a^2-2*mu*ln(a*y*z)*a^2+la*ln(a*y*z)^2*a^2-4*a*mu*x+4*a*x*la*ln(a*y*z)-3*a^2*la*ln(a*y*z)+x^2*mu*a^2+mu*x^2+x^2*la-x^2*la*ln(a*y*z)-2*x*a*la+a^2*la)/a^2;
    }
    else if ((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        return 1/2*(x^2*mu*a^2+mu*z^2*a^2-2*mu*ln(x*a*z)*a^2+la*ln(x*a*z)^2*a^2-4*y*mu*a+4*a*y*la*ln(x*a*z)-3*a^2*la*ln(x*a*z)+mu*y^2*a^2+mu*y^2+y^2*la-y^2*la*ln(x*a*z)-2*y*a*la+a^2*la)/a^2;
    }
    else if ((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        return 1/2*(x^2*mu*a^2+mu*y^2*a^2-2*mu*ln(x*y*a)*a^2+la*ln(x*y*a)^2*a^2-4*z*mu*a+4*a*z*la*ln(x*y*a)-3*a^2*la*ln(x*y*a)+mu*z^2*a^2+mu*z^2+z^2*la-z^2*la*ln(x*y*a)-2*z*a*la+a^2*la)/a^2;
    }
    else if ((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        return 1/4*(8*a^3*y*la*ln(z*a^2)-2*a^2*y^2*la*ln(z*a^2)+8*x*a^3*la*ln(z*a^2)+16*x*y*a^2*la-2*x^2*a^2*la*ln(z*a^2)+2*mu*y^2*a^2-8*x*mu*a^3+2*x^2*mu*a^2+6*mu*a^4-4*y^2*x*a*la+5*y^2*a^2*la+y^2*x^2*la-8*mu*y*a^3-4*x^2*y*a*la+5*x^2*a^2*la+13*a^4*la+2*a^4*mu*z^2-4*a^4*mu*ln(z*a^2)+2*a^4*la*ln(z*a^2)^2-12*a^4*la*ln(z*a^2)+2*a^4*mu*y^2-16*a^3*y*la-16*x*a^3*la+2*x^2*mu*a^4)/a^4;
    }
    else if ((dx < 0) && (dy >= 0) && (dz < 0)) // Rxz
    {
        return 1/4*(2*mu*z^2*a^2-8*x*mu*a^3+2*x^2*mu*a^2+6*mu*a^4-4*z^2*x*a*la+5*z^2*a^2*la+z^2*x^2*la+5*x^2*a^2*la-8*mu*z*a^3-4*x^2*z*a*la+13*a^4*la+2*a^4*mu*z^2+2*a^4*mu*y^2-16*x*a^3*la+2*x^2*mu*a^4-16*z*a^3*la+16*z*x*a^2*la-4*a^4*mu*ln(y*a^2)+2*a^4*la*ln(y*a^2)^2-12*a^4*la*ln(y*a^2)+8*a^3*z*la*ln(y*a^2)-2*a^2*z^2*la*ln(y*a^2)+8*x*a^3*la*ln(y*a^2)-2*x^2*a^2*la*ln(y*a^2))/a^4;
    }
    else if ((dx >= 0) && (dy < 0) && (dz < 0)) // Rzy
    {
        return 1/4*(2*mu*y^2*a^2+2*mu*z^2*a^2+6*mu*a^4+5*y^2*a^2*la+5*z^2*a^2*la-8*mu*y*a^3-4*z^2*y*a*la+z^2*y^2*la-8*mu*z*a^3-4*y^2*z*a*la+13*a^4*la+2*a^4*mu*z^2+2*a^4*mu*y^2-16*a^3*y*la+2*x^2*mu*a^4-16*z*a^3*la+16*z*y*a^2*la-4*a^4*mu*ln(x*a^2)+2*a^4*la*ln(x*a^2)^2-12*a^4*la*ln(x*a^2)+8*a^3*y*la*ln(x*a^2)-2*a^2*y^2*la*ln(x*a^2)+8*z*a^3*la*ln(x*a^2)-2*z^2*a^2*la*ln(x*a^2))/a^4;
    }
    else // Rxyz
    {
        return 1/4*(16*x*y*a^2*la+2*mu*y^2*a^2+2*mu*z^2*a^2-8*x*mu*a^3+2*x^2*mu*a^2+12*mu*a^4-4*y^2*x*a*la-4*z^2*x*a*la+8*y^2*a^2*la+y^2*x^2*la+8*z^2*a^2*la+z^2*x^2*la-8*mu*y*a^3-4*x^2*y*a*la-4*z^2*y*a*la+8*x^2*a^2*la+z^2*y^2*la-8*mu*z*a^3-4*x^2*z*a*la-4*y^2*z*a*la+33*a^4*la+2*a^4*mu*z^2+2*a^4*mu*y^2-28*a^3*y*la-28*x*a^3*la+2*x^2*mu*a^4-28*z*a^3*la+16*z*y*a^2*la+16*z*x*a^2*la-4*a^4*mu*ln(a^3)+2*a^4*la*ln(a^3)^2-18*a^4*la*ln(a^3)+8*a^3*y*la*ln(a^3)-2*a^2*y^2*la*ln(a^3)+8*z*a^3*la*ln(a^3)-2*z^2*a^2*la*ln(a^3)+8*x*a^3*la*ln(a^3)-2*x^2*a^2*la*ln(a^3))/a^4;
    }
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

    T mu = constant_mu;
    T la = constant_lambda;
    
    T a = extrapolation_cutoff;
    
    if ((dx >= 0) && (dy >= 0))
    {  
        T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda,J=F.Determinant();
        return scale_mu*F-(scale_mu-scale_lambda*log(J))*F.Inverse();
    }
    else if ((dx < 0) && (dy >= 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = (-2*mu*a+2*a*la*log(y*a)-a*la+x*mu*sqr(a)+mu*x+la*x-la*log(y*a)*x)/sqr(a);
        result.x22 = 0.5*(2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(y*a)-3*la*sqr(a)+4*x*a*la-sqr(x)*la)/y/sqr(a);
        return scale*result;
    }
    else if ((dx >= 0) && (dy < 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = 0.5*(2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*a)-3*la*sqr(a)+4*y*a*la-sqr(y)*la)/x/sqr(a);
        result.x22 = (-2*mu*a+2*a*la*log(x*a)-a*la+mu*sqr(a)*y+mu*y+la*y-la*log(x*a)*y)/sqr(a);
        return scale*result;
    }
    else // ((dx < 0) && (dy < 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = 0.5*(-4*mu*cube(a)+2*x*mu*sqr(a)+4*cube(a)*la*log(sqr(a))+8*y*la*sqr(a)-2*sqr(y)*a*la-2*x*sqr(a)*la*log(sqr(a))-4*x*y*a*la-8*la*cube(a)+5*x*la*sqr(a)+2*x*mu*(a*a*a*a)+sqr(y)*la*x)/(a*a*a*a);
        result.x22 = 0.5*(2*mu*sqr(a)*y-4*mu*cube(a)+4*cube(a)*la*log(sqr(a))-2*sqr(a)*la*log(sqr(a))*y+8*x*la*sqr(a)-4*x*y*a*la-2*sqr(x)*a*la+5*y*la*sqr(a)+2*mu*(a*a*a*a)*y-8*la*cube(a)+sqr(x)*la*y)/(a*a*a*a);
        return scale*result;
    }
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const
{
    T x = F.x11;
    T y = F.x22;
    T z = F.x22;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = constant_mu;
    T la = constant_lambda;
    
    T a = extrapolation_cutoff;
     
    if ((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda,J=F.Determinant();
        return scale_mu*F-(scale_mu-scale_lambda*log(J))*F.Inverse();
    }
    else if ((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=(-2*mu*a+2*a*la*ln(a*y*z)+x*mu*a^2+mu*x+x*la-x*la*ln(a*y*z)-a*la)/a^2;
        result.x22=1/2*(2*mu*y^2*a^2-2*mu*a^2+2*a^2*la*ln(a*y*z)+4*x*a*la-3*a^2*la-x^2*la)/y/a^2;
        result.x33=1/2*(2*mu*z^2*a^2-2*mu*a^2+2*a^2*la*ln(a*y*z)+4*x*a*la-3*a^2*la-x^2*la)/z/a^2;
        return scale*result;
    }
    else if ((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=1/2*(2*x^2*mu*a^2-2*mu*a^2+2*a^2*la*ln(x*a*z)+4*y*a*la-3*a^2*la-y^2*la)/x/a^2;
        result.x22=(-2*mu*a+2*a*la*ln(x*a*z)+mu*y*a^2+mu*y+y*la-y*la*ln(x*a*z)-a*la)/a^2;
        result.x33=1/2*(2*mu*z^2*a^2-2*mu*a^2+2*a^2*la*ln(x*a*z)+4*y*a*la-3*a^2*la-y^2*la)/z/a^2;
        return scale*result;
    }
    else if ((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=1/2*(2*x^2*mu*a^2-2*mu*a^2+2*a^2*la*ln(x*y*a)+4*z*a*la-3*a^2*la-z^2*la)/x/a^2;
        result.x22=1/2*(2*mu*y^2*a^2-2*mu*a^2+2*a^2*la*ln(x*y*a)+4*z*a*la-3*a^2*la-z^2*la)/y/a^2;
        result.x33=(-2*mu*a+2*a*la*ln(x*y*a)+mu*z*a^2+mu*z+z*la-z*la*ln(x*y*a)-a*la)/a^2;
        return scale*result;
    }
    else if ((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=1/2*(4*a^3*la*ln(z*a^2)+8*y*a^2*la-2*x*a^2*la*ln(z*a^2)-4*mu*a^3+2*x*mu*a^2-2*y^2*a*la+y^2*x*la-4*x*y*a*la+5*x*a^2*la-8*a^3*la+2*x*mu*a^4)/a^4;
        result.x22=1/2*(4*a^3*la*ln(z*a^2)-2*a^2*y*la*ln(z*a^2)+8*x*a^2*la+2*mu*y*a^2-4*x*y*a*la+5*y*a^2*la+x^2*y*la-4*mu*a^3-2*x^2*a*la+2*a^4*mu*y-8*a^3*la)/a^4;
        result.x33=1/2/a^2*(4*y*a*la-y^2*la+4*x*a*la-x^2*la+2*mu*z^2*a^2-2*mu*a^2+2*a^2*la*ln(z*a^2)-6*a^2*la)/z;
        return scale*result;
    }
    else if ((dx < 0) && (dy >= 0) && (dz < 0)) // Rxz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=1/2*(-4*mu*a^3+2*x*mu*a^2-2*z^2*a*la+z^2*x*la+5*x*a^2*la-4*x*z*a*la-8*a^3*la+2*x*mu*a^4+8*z*a^2*la+4*a^3*la*ln(y*a^2)-2*x*a^2*la*ln(y*a^2))/a^4;
        result.x22=1/2/a^2*(2*mu*y^2*a^2-2*mu*a^2+2*a^2*la*ln(y*a^2)-6*a^2*la+4*z*a*la-z^2*la+4*x*a*la-x^2*la)/y;
        result.x33=1/2*(2*mu*z*a^2-4*x*z*a*la+5*z*a^2*la+x^2*z*la-4*mu*a^3-2*x^2*a*la+2*a^4*mu*z-8*a^3*la+8*x*a^2*la+4*a^3*la*ln(y*a^2)-2*a^2*z*la*ln(y*a^2))/a^4;
        return scale*result;
    }
    else if ((dx >= 0) && (dy < 0) && (dz < 0)) // Rzy
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=1/2/a^2*(2*x^2*mu*a^2-2*mu*a^2+2*a^2*la*ln(x*a^2)-6*a^2*la+4*y*a*la-y^2*la+4*z*a*la-z^2*la)/x;
        result.x22=1/2*(2*mu*y*a^2+5*y*a^2*la-4*mu*a^3-2*z^2*a*la+z^2*y*la-4*z*y*a*la+2*a^4*mu*y-8*a^3*la+8*z*a^2*la+4*a^3*la*ln(x*a^2)-2*a^2*y*la*ln(x*a^2))/a^4;
        result.x33=1/2*(2*mu*z*a^2+5*z*a^2*la-4*z*y*a*la+y^2*z*la-4*mu*a^3-2*y^2*a*la+2*a^4*mu*z-8*a^3*la+8*y*a^2*la+4*a^3*la*ln(x*a^2)-2*z*a^2*la*ln(x*a^2))/a^4;
        return scale*result;
    }
    else // Rxyz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=1/2*(8*y*a^2*la-4*mu*a^3+2*x*mu*a^2-2*y^2*a*la-2*z^2*a*la+y^2*x*la+z^2*x*la-4*x*y*a*la+8*x*a^2*la-4*x*z*a*la-14*a^3*la+2*x*mu*a^4+8*z*a^2*la+4*a^3*la*ln(a^3)-2*x*a^2*la*ln(a^3))/a^4;
        result.x22=1/2*(8*x*a^2*la+2*mu*y*a^2-4*x*y*a*la+8*y*a^2*la+x^2*y*la-4*mu*a^3-2*x^2*a*la-2*z^2*a*la+z^2*y*la-4*z*y*a*la+2*a^4*mu*y-14*a^3*la+8*z*a^2*la+4*a^3*la*ln(a^3)-2*a^2*y*la*ln(a^3))/a^4;
        result.x33=1/2*(2*mu*z*a^2-4*x*z*a*la+8*z*a^2*la+x^2*z*la-4*z*y*a*la+y^2*z*la-4*mu*a^3-2*x^2*a*la-2*y^2*a*la+2*a^4*mu*z-14*a^3*la+8*y*a^2*la+8*x*a^2*la+4*a^3*la*ln(a^3)-2*z*a^2*la*ln(a^3))/a^4;
        return scale*result;
    }
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

    T mu = constant_mu;
    T la = constant_lambda;
    
    T a = extrapolation_cutoff;
    
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
        dP_dF.x1111=-(-mu*sqr(a)-mu-la+la*log(y*a))/sqr(a);
        dP_dF.x2222=-0.5*(-2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)-5*la*sqr(a)+2*sqr(a)*la*log(y*a)+4*x*a*la-sqr(x)*la)/sqr(y)/sqr(a);
        dP_dF.x2211=la*(2*a-x)/y/sqr(a);
  
        T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T xmy = x-y; if (fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;

        dP_dF.x2112=-0.5*(-2*x*mu*sqr(a)+2*x*sqr(a)*la*log(y*a)-3*x*la*sqr(a)+4*sqr(x)*a*la-cube(x)*la+4*sqr(y)*mu*a-4*sqr(y)*a*la*log(y*a)+2*sqr(y)*a*la-2*sqr(y)*mu*x-2*sqr(y)*la*x+2*sqr(y)*la*log(y*a)*x)/y/sqr(a)/(sqr(x)-sqr(y));
        dP_dF.x2121=-0.5*(2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(y*a)-3*la*sqr(a)+6*x*a*la-3*sqr(x)*la+4*mu*x*a-4*a*la*log(y*a)*x-2*sqr(x)*mu*sqr(a)-2*mu*sqr(x)+2*la*log(y*a)*sqr(x))/sqr(a)/(sqr(x)-sqr(y));

    }
    else if ((dx >= 0) && (dy < 0))
    {
        dP_dF.x1111=-0.5*(-2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)-5*la*sqr(a)+2*sqr(a)*la*log(x*a)+4*y*a*la-sqr(y)*la)/sqr(x)/sqr(a);
        dP_dF.x2222=-(-mu*sqr(a)-mu-la+la*log(x*a))/sqr(a);
        dP_dF.x2211=la*(2*a-y)/x/sqr(a);
 
        T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T xmy = x-y; if (fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;

        dP_dF.x2112=0.5*(4*sqr(x)*mu*a-4*sqr(x)*a*la*log(x*a)+2*sqr(x)*a*la-2*sqr(x)*mu*y-2*sqr(x)*la*y+2*sqr(x)*la*log(x*a)*y-2*mu*sqr(a)*y+2*y*sqr(a)*la*log(x*a)-3*y*la*sqr(a)+4*sqr(y)*a*la-cube(y)*la)/x/sqr(a)/(sqr(x)-sqr(y));
        dP_dF.x2121=0.5*(4*mu*y*a-4*a*y*la*log(x*a)+6*y*a*la-2*mu*sqr(y)*sqr(a)-2*mu*sqr(y)-3*sqr(y)*la+2*sqr(y)*la*log(x*a)+2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*a)-3*la*sqr(a))/sqr(a)/(sqr(x)-sqr(y));
    }
    else // ((dx < 0) && (dy < 0))
    {
        dP_dF.x1111=-0.5*(-5*la*sqr(a)+2*la*sqr(a)*log(sqr(a))-2*mu*(a*a*a*a)-2*mu*sqr(a)-sqr(y)*la+4*y*a*la)/(a*a*a*a);
        dP_dF.x2222=-0.5*(-2*mu*sqr(a)+2*la*sqr(a)*log(sqr(a))+4*x*a*la-5*la*sqr(a)-2*mu*(a*a*a*a)-sqr(x)*la)/(a*a*a*a);
        dP_dF.x2211=la*(4*sqr(a)-2*y*a-2*x*a+x*y)/(a*a*a*a);

        T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
  
        dP_dF.x2112=-0.5*(sqr(x)*la*y-2*sqr(x)*a*la+8*x*la*sqr(a)-6*x*y*a*la+sqr(y)*la*x+4*cube(a)*la*log(sqr(a))-2*sqr(y)*a*la-4*mu*cube(a)-8*la*cube(a)+8*y*la*sqr(a))/(x+y)/(a*a*a*a);
        dP_dF.x2121=0.5*(2*mu*x*a-2*x*a*la*log(sqr(a))-2*y*la*x+5*x*a*la+2*x*mu*cube(a)-4*mu*sqr(a)-8*la*sqr(a)+4*la*sqr(a)*log(sqr(a))+2*mu*y*a-2*a*y*la*log(sqr(a))+5*y*a*la+2*mu*y*cube(a))/(x+y)/cube(a);
    }

    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int triangle) const
{
    T x = F.x11;
    T y = F.x22;
    T z = F.x22;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = constant_mu;
    T la = constant_lambda;
    
    T a = extrapolation_cutoff;
     
    if ((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        DIAGONAL_MATRIX<T,3> F_inverse=F.Inverse();
        T mu_minus_lambda_logJ=constant_mu+constant_lambda*log(F_inverse.Determinant());
        SYMMETRIC_MATRIX<T,3> F_inverse_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(F_inverse.To_Vector());
        dPi_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;
        dPi_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
        dPi_dF.x3333=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x33;
        dPi_dF.x2211=constant_lambda*F_inverse_outer.x21;
        dPi_dF.x3311=constant_lambda*F_inverse_outer.x31;
        dPi_dF.x3322=constant_lambda*F_inverse_outer.x32;
        dPi_dF.x2121=constant_mu;
        dPi_dF.x3131=constant_mu;
        dPi_dF.x3232=constant_mu;
        dPi_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;
        dPi_dF.x3113=mu_minus_lambda_logJ*F_inverse_outer.x31;
        dPi_dF.x3223=mu_minus_lambda_logJ*F_inverse_outer.x32;
    }
    else if ((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
        dPi_dF.x1111=-(-mu*a^2-mu-la+la*ln(a*y*z))/a^2;
        dPi_dF.x2222=-1/2*(-2*mu*y^2*a^2-2*mu*a^2-5*a^2*la+2*a^2*la*ln(a*y*z)+4*x*a*la-x^2*la)/y^2/a^2;
        dPi_dF.x3333=-1/2*(-2*mu*z^2*a^2-2*mu*a^2-5*a^2*la+2*a^2*la*ln(a*y*z)+4*x*a*la-x^2*la)/z^2/a^2;
        dPi_dF.x2211=la*(2*a-x)/y/a^2;
        dPi_dF.x3311=la*(2*a-x)/z/a^2;
        dPi_dF.x3322=la/y/z;
        dPi_dF.x2121=-1/2*(2*mu*y^2*a^2-2*mu*a^2+2*a^2*la*ln(a*y*z)+6*x*a*la-3*a^2*la-3*x^2*la+4*a*mu*x-4*a*x*la*ln(a*y*z)-2*x^2*mu*a^2-2*mu*x^2+2*x^2*la*ln(a*y*z))/a^2/(x^2-y^2);
        dPi_dF.x3232=mu;
        dPi_dF.x2112=-1/2*(-2*x*mu*a^2+2*x*a^2*la*ln(a*y*z)+4*x^2*a*la-3*x*a^2*la-x^3*la+4*y^2*mu*a-4*y^2*a*la*ln(a*y*z)-2*y^2*mu*x-2*y^2*x*la+2*y^2*x*la*ln(a*y*z)+2*y^2*a*la)/y/a^2/(x^2-y^2);
        dPi_dF.x3113=-1/2*(-2*x*mu*a^2+2*x*a^2*la*ln(a*y*z)+4*x^2*a*la-3*x*a^2*la-x^3*la+4*z^2*mu*a-4*z^2*a*la*ln(a*y*z)-2*z^2*mu*x-2*z^2*x*la+2*z^2*x*la*ln(a*y*z)+2*z^2*a*la)/z/a^2/(x^2-z^2);
        dPi_dF.x3223=-1/2*(-2*mu*a^2+4*x*a*la-3*a^2*la-x^2*la+2*a^2*la*ln(a*y*z))/z/a^2/y;
        dPi_dF.x3131=-1/2*(2*mu*z^2*a^2-2*mu*a^2+2*a^2*la*ln(a*y*z)+6*x*a*la-3*a^2*la-3*x^2*la+4*a*mu*x-4*a*x*la*ln(a*y*z)-2*x^2*mu*a^2-2*mu*x^2+2*x^2*la*ln(a*y*z))/a^2/(x^2-z^2);
    }
    else if ((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        dPi_dF.x1111=-1/2*(-2*x^2*mu*a^2-2*mu*a^2-5*a^2*la+2*a^2*la*ln(x*a*z)+4*y*a*la-y^2*la)/x^2/a^2;
        dPi_dF.x2222=-(-mu*a^2-mu-la+la*ln(x*a*z))/a^2;
        dPi_dF.x3333=-1/2*(-2*mu*z^2*a^2-2*mu*a^2-5*a^2*la+2*a^2*la*ln(x*a*z)+4*y*a*la-y^2*la)/z^2/a^2;
        dPi_dF.x2211=la*(2*a-y)/x/a^2;
        dPi_dF.x3311=la/x/z;
        dPi_dF.x3322=la*(2*a-y)/z/a^2;
        dPi_dF.x2121=1/2*(4*y*mu*a-4*a*y*la*ln(x*a*z)-2*mu*y^2*a^2-2*mu*y^2-3*y^2*la+2*y^2*la*ln(x*a*z)+6*y*a*la+2*x^2*mu*a^2-2*mu*a^2+2*a^2*la*ln(x*a*z)-3*a^2*la)/a^2/(x^2-y^2);
        dPi_dF.x3232=-1/2*(2*mu*z^2*a^2-2*mu*a^2+2*a^2*la*ln(x*a*z)+6*y*a*la-3*a^2*la-3*y^2*la+4*y*mu*a-4*a*y*la*ln(x*a*z)-2*mu*y^2*a^2-2*mu*y^2+2*y^2*la*ln(x*a*z))/a^2/(y^2-z^2);
        dPi_dF.x2112=1/2*(4*x^2*mu*a-4*x^2*a*la*ln(x*a*z)-2*x^2*mu*y-2*x^2*y*la+2*x^2*y*la*ln(x*a*z)+2*x^2*a*la-2*mu*y*a^2+2*y*a^2*la*ln(x*a*z)+4*y^2*a*la-3*y*a^2*la-y^3*la)/x/a^2/(x^2-y^2);
        dPi_dF.x3113=-1/2*(-2*mu*a^2+4*y*a*la-3*a^2*la-y^2*la+2*a^2*la*ln(x*a*z))/z/a^2/x;
        dPi_dF.x3223=-1/2*(-2*mu*y*a^2+2*y*a^2*la*ln(x*a*z)+4*y^2*a*la-3*y*a^2*la-y^3*la+4*z^2*mu*a-4*z^2*a*la*ln(x*a*z)-2*z^2*mu*y-2*z^2*y*la+2*z^2*y*la*ln(x*a*z)+2*z^2*a*la)/z/a^2/(y^2-z^2);
        dPi_dF.x3131=mu;
    }
    else if ((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        dPi_dF.x1111=-1/2*(-2*x^2*mu*a^2-2*mu*a^2-5*a^2*la+2*a^2*la*ln(x*y*a)+4*z*a*la-z^2*la)/x^2/a^2;
        dPi_dF.x2222=-1/2*(-2*mu*y^2*a^2-2*mu*a^2-5*a^2*la+2*a^2*la*ln(x*y*a)+4*z*a*la-z^2*la)/y^2/a^2;
        dPi_dF.x3333=-(-mu*a^2-mu-la+la*ln(x*y*a))/a^2;
        dPi_dF.x2211=la/y/x;
        dPi_dF.x3311=la*(2*a-z)/x/a^2;
        dPi_dF.x3322=la*(2*a-z)/y/a^2;
        dPi_dF.x2121=mu;
        dPi_dF.x3232=1/2*(4*z*mu*a-4*a*z*la*ln(x*y*a)-2*mu*z^2*a^2-2*mu*z^2-3*z^2*la+2*z^2*la*ln(x*y*a)+6*z*a*la+2*mu*y^2*a^2-2*mu*a^2+2*a^2*la*ln(x*y*a)-3*a^2*la)/a^2/(y^2-z^2);
        dPi_dF.x2112=-1/2*(-2*mu*a^2+4*z*a*la-3*a^2*la-z^2*la+2*a^2*la*ln(x*y*a))/y/a^2/x;
        dPi_dF.x3113=1/2*(4*x^2*mu*a-4*x^2*a*la*ln(x*y*a)-2*x^2*mu*z-2*x^2*z*la+2*x^2*z*la*ln(x*y*a)+2*x^2*a*la-2*mu*z*a^2+2*z*a^2*la*ln(x*y*a)+4*z^2*a*la-3*z*a^2*la-z^3*la)/x/a^2/(x^2-z^2);
        dPi_dF.x3223=1/2*(4*y^2*mu*a-4*y^2*a*la*ln(x*y*a)-2*y^2*mu*z-2*y^2*z*la+2*y^2*z*la*ln(x*y*a)+2*y^2*a*la-2*mu*z*a^2+2*z*a^2*la*ln(x*y*a)+4*z^2*a*la-3*z*a^2*la-z^3*la)/y/a^2/(y^2-z^2);
        dPi_dF.x3131=1/2*(4*z*mu*a-4*a*z*la*ln(x*y*a)-2*mu*z^2*a^2-2*mu*z^2-3*z^2*la+2*z^2*la*ln(x*y*a)+6*z*a*la+2*x^2*mu*a^2-2*mu*a^2+2*a^2*la*ln(x*y*a)-3*a^2*la)/a^2/(x^2-z^2);
    }
    else if ((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        dPi_dF.x1111=-1/2*(-5*a^2*la+2*a^2*la*ln(z*a^2)-2*mu*a^4-2*mu*a^2-y^2*la+4*y*a*la)/a^4;
        dPi_dF.x2222=-1/2*(2*a^2*la*ln(z*a^2)-2*mu*a^2+4*x*a*la-5*a^2*la-x^2*la-2*mu*a^4)/a^4;
        dPi_dF.x3333=-1/2/a^2*(4*y*a*la-y^2*la+4*x*a*la-x^2*la-2*mu*z^2*a^2-2*mu*a^2-8*a^2*la+2*a^2*la*ln(z*a^2))/z^2;
        dPi_dF.x2211=la*(4*a^2-2*a*y+x*y-2*x*a)/a^4;
        dPi_dF.x3311=la*(2*a-x)/z/a^2;
        dPi_dF.x3322=la*(2*a-y)/z/a^2;
        dPi_dF.x2121=1/2*(-2*y*x*la+2*a*mu*x+5*x*a*la-2*x*a*la*ln(z*a^2)+2*x*mu*a^3-8*a^2*la+4*a^2*la*ln(z*a^2)-4*mu*a^2+2*y*mu*a+5*y*a*la-2*a*y*la*ln(z*a^2)+2*mu*y*a^3)/(x+y)/a^3;
        dPi_dF.x3232=-1/2*(12*a^3*y*la-6*y^2*a^2*la+4*x*a^3*la-x^2*a^2*la+2*a^4*mu*z^2-2*mu*a^4+2*a^4*la*ln(z*a^2)-6*a^4*la-4*a^3*y*la*ln(z*a^2)+2*a^2*y^2*la*ln(z*a^2)-8*x*y*a^2*la-2*mu*y^2*a^2+4*y^2*x*a*la-y^2*x^2*la+4*mu*y*a^3+2*x^2*y*a*la-2*a^4*mu*y^2)/a^4/(y^2-z^2);
        dPi_dF.x2112=-1/2*(x^2*y*la-2*x^2*a*la+8*x*a^2*la-6*x*y*a*la+y^2*x*la+4*a^3*la*ln(z*a^2)-8*a^3*la-4*mu*a^3-2*y^2*a*la+8*y*a^2*la)/(x+y)/a^4;
        dPi_dF.x3113=-1/2*(4*x*a^3*y*la-x*a^2*y^2*la+4*x^2*a^3*la-x^3*a^2*la-2*x*mu*a^4+2*x*a^4*la*ln(z*a^2)-6*x*a^4*la-4*z^2*a^3*la*ln(z*a^2)-8*z^2*y*a^2*la+2*z^2*x*a^2*la*ln(z*a^2)+4*z^2*mu*a^3-2*x*mu*z^2*a^2+2*z^2*y^2*a*la-z^2*y^2*x*la+4*z^2*x*y*a*la-5*z^2*x*a^2*la+8*z^2*a^3*la)/z/a^4/(x^2-z^2);
        dPi_dF.x3223=-1/2*(4*y^2*a^3*la-y^3*a^2*la+4*x*a^3*y*la-y*a^2*x^2*la-2*a^4*mu*y+2*y*a^4*la*ln(z*a^2)-6*y*a^4*la-4*z^2*a^3*la*ln(z*a^2)+2*z^2*a^2*y*la*ln(z*a^2)-8*z^2*x*a^2*la-2*y*mu*z^2*a^2+4*z^2*x*y*a*la-5*z^2*y*a^2*la-z^2*x^2*y*la+4*z^2*mu*a^3+2*z^2*x^2*a*la+8*z^2*a^3*la)/z/a^4/(y^2-z^2);
        dPi_dF.x3131=-1/2*(4*a^3*y*la-y^2*a^2*la+12*x*a^3*la-6*x^2*a^2*la+2*a^4*mu*z^2-2*mu*a^4+2*a^4*la*ln(z*a^2)-6*a^4*la-4*x*a^3*la*ln(z*a^2)-8*x*y*a^2*la+2*x^2*a^2*la*ln(z*a^2)+4*x*mu*a^3-2*x^2*mu*a^2+2*y^2*x*a*la-y^2*x^2*la+4*x^2*y*a*la-2*x^2*mu*a^4)/a^4/(x^2-z^2);
    }
    else if ((dx < 0) && (dy >= 0) && (dz < 0)) // Rxz
    {
        dPi_dF.x1111=-1/2*(-5*a^2*la+2*a^2*la*ln(y*a^2)-2*mu*a^4-2*mu*a^2-z^2*la+4*z*a*la)/a^4;
        dPi_dF.x2222=-1/2/a^2*(-2*mu*y^2*a^2-2*mu*a^2-8*a^2*la+2*a^2*la*ln(y*a^2)+4*z*a*la-z^2*la+4*x*a*la-x^2*la)/y^2;
        dPi_dF.x3333=-1/2*(-2*mu*a^2+4*x*a*la-5*a^2*la-x^2*la-2*mu*a^4+2*a^2*la*ln(y*a^2))/a^4;
        dPi_dF.x2211=la*(2*a-x)/y/a^2;
        dPi_dF.x3311=la*(-2*a*z+x*z-2*x*a+4*a^2)/a^4;
        dPi_dF.x3322=la*(2*a-z)/y/a^2;
        dPi_dF.x2121=-1/2*(2*a^4*mu*y^2-2*mu*a^4+2*a^4*la*ln(y*a^2)-6*a^4*la+4*z*a^3*la-z^2*a^2*la+12*x*a^3*la-6*x^2*a^2*la+4*x*mu*a^3-2*x^2*mu*a^2+2*z^2*x*a*la-z^2*x^2*la+4*x^2*z*a*la-2*x^2*mu*a^4-8*z*x*a^2*la-4*x*a^3*la*ln(y*a^2)+2*x^2*a^2*la*ln(y*a^2))/a^4/(x^2-y^2);
        dPi_dF.x3232=1/2*(-2*mu*z^2*a^2+4*z^2*x*a*la-6*z^2*a^2*la-z^2*x^2*la+4*mu*z*a^3+2*x^2*z*a*la-2*a^4*mu*z^2+12*z*a^3*la-8*z*x*a^2*la-4*a^3*z*la*ln(y*a^2)+2*a^2*z^2*la*ln(y*a^2)+2*a^4*mu*y^2-2*mu*a^4+2*a^4*la*ln(y*a^2)-6*a^4*la+4*x*a^3*la-x^2*a^2*la)/a^4/(y^2-z^2);
        dPi_dF.x2112=-1/2*(-2*x*mu*a^4+2*x*a^4*la*ln(y*a^2)-6*x*a^4*la+4*x*a^3*z*la-z^2*x*a^2*la+4*x^2*a^3*la-x^3*a^2*la+4*a^3*mu*y^2-2*x*mu*y^2*a^2+2*z^2*y^2*a*la-z^2*y^2*x*la-5*x*a^2*y^2*la+4*y^2*x*z*a*la+8*y^2*a^3*la-8*y^2*z*a^2*la-4*y^2*a^3*la*ln(y*a^2)+2*y^2*x*a^2*la*ln(y*a^2))/y/a^4/(x^2-y^2);
        dPi_dF.x3113=-1/2*(-2*x^2*a*la+x^2*z*la+8*x*a^2*la-6*x*z*a*la+z^2*x*la-4*mu*a^3+4*a^3*la*ln(y*a^2)-2*z^2*a*la-8*a^3*la+8*z*a^2*la)/(x+z)/a^4;
        dPi_dF.x3223=1/2*(-2*y^2*mu*z*a^2+4*y^2*x*z*a*la-5*y^2*z*a^2*la-x^2*y^2*z*la+4*a^3*mu*y^2+2*x^2*y^2*a*la+8*y^2*a^3*la-8*x*a^2*y^2*la-4*y^2*a^3*la*ln(y*a^2)+2*y^2*a^2*z*la*ln(y*a^2)-2*a^4*mu*z+2*z*a^4*la*ln(y*a^2)-6*z*a^4*la+4*z^2*a^3*la-z^3*a^2*la+4*x*a^3*z*la-x^2*z*a^2*la)/y/a^4/(y^2-z^2);
        dPi_dF.x3131=1/2*(5*x*a*la-2*z*x*la+2*x*mu*a^3+2*a*mu*x-2*x*a*la*ln(y*a^2)-8*a^2*la-4*mu*a^2+4*a^2*la*ln(y*a^2)+5*z*a*la+2*mu*z*a^3+2*z*mu*a-2*a*z*la*ln(y*a^2))/(x+z)/a^3;
    }
    else if ((dx >= 0) && (dy < 0) && (dz < 0)) // Rzy
    {
        dPi_dF.x1111=-1/2/a^2*(-2*x^2*mu*a^2-2*mu*a^2-8*a^2*la+2*a^2*la*ln(x*a^2)+4*y*a*la-y^2*la+4*z*a*la-z^2*la)/x^2;
        dPi_dF.x2222=-1/2*(-2*mu*a^2-5*a^2*la-z^2*la+4*z*a*la-2*mu*a^4+2*a^2*la*ln(x*a^2))/a^4;
        dPi_dF.x3333=-1/2*(-5*a^2*la+2*a^2*la*ln(x*a^2)-2*mu*a^4-2*mu*a^2-y^2*la+4*y*a*la)/a^4;
        dPi_dF.x2211=la*(2*a-y)/x/a^2;
        dPi_dF.x3311=la*(2*a-z)/x/a^2;
        dPi_dF.x3322=la*(-2*a*z+y*z-2*a*y+4*a^2)/a^4;
        dPi_dF.x2121=1/2*(-2*mu*y^2*a^2-6*y^2*a^2*la+4*mu*y*a^3+2*z^2*y*a*la-z^2*y^2*la+4*y^2*z*a*la-2*a^4*mu*y^2+12*a^3*y*la-8*z*y*a^2*la-4*a^3*y*la*ln(x*a^2)+2*a^2*y^2*la*ln(x*a^2)+2*x^2*mu*a^4-2*mu*a^4+2*a^4*la*ln(x*a^2)-6*a^4*la+4*z*a^3*la-z^2*a^2*la)/a^4/(x^2-y^2);
        dPi_dF.x3232=1/2*(-2*z*y*la+2*mu*y*a^3+2*y*mu*a+5*y*a*la-2*a*y*la*ln(x*a^2)-8*a^2*la+4*a^2*la*ln(x*a^2)-4*mu*a^2+2*mu*z*a^3+2*z*mu*a+5*z*a*la-2*z*a*la*ln(x*a^2))/(y+z)/a^3;
        dPi_dF.x2112=1/2*(-2*x^2*mu*y*a^2-5*y*a^2*x^2*la+4*x^2*mu*a^3+2*z^2*x^2*a*la-z^2*x^2*y*la+4*x^2*z*y*a*la+8*x^2*a^3*la-8*x^2*z*a^2*la-4*x^2*a^3*la*ln(x*a^2)+2*x^2*a^2*y*la*ln(x*a^2)-2*a^4*mu*y+2*y*a^4*la*ln(x*a^2)-6*y*a^4*la+4*y^2*a^3*la-y^3*a^2*la+4*y*a^3*z*la-z^2*y*a^2*la)/x/a^4/(x^2-y^2);
        dPi_dF.x3113=1/2*(-2*x^2*mu*z*a^2-5*x^2*z*a^2*la+4*x^2*z*y*a*la-x^2*y^2*z*la+4*x^2*mu*a^3+2*x^2*y^2*a*la+8*x^2*a^3*la-8*y*a^2*x^2*la-4*x^2*a^3*la*ln(x*a^2)+2*x^2*z*a^2*la*ln(x*a^2)-2*a^4*mu*z+2*z*a^4*la*ln(x*a^2)-6*z*a^4*la+4*y*a^3*z*la-y^2*z*a^2*la+4*z^2*a^3*la-z^3*a^2*la)/x/a^4/(x^2-z^2);
        dPi_dF.x3223=-1/2*(-2*y^2*a*la+y^2*z*la+8*y*a^2*la-6*z*y*a*la+z^2*y*la-4*mu*a^3+4*a^3*la*ln(x*a^2)-2*z^2*a*la-8*a^3*la+8*z*a^2*la)/(y+z)/a^4;
        dPi_dF.x3131=1/2*(-2*mu*z^2*a^2-6*z^2*a^2*la+4*z^2*y*a*la-z^2*y^2*la+4*mu*z*a^3+2*y^2*z*a*la-2*a^4*mu*z^2+12*z*a^3*la-8*z*y*a^2*la-4*z*a^3*la*ln(x*a^2)+2*z^2*a^2*la*ln(x*a^2)+2*x^2*mu*a^4-2*mu*a^4+2*a^4*la*ln(x*a^2)-6*a^4*la+4*a^3*y*la-y^2*a^2*la)/a^4/(x^2-z^2);
    }
    else // Rxyz
    {
        dPi_dF.x1111=-1/2*(-8*a^2*la+2*a^2*la*ln(a^3)-2*mu*a^4-2*mu*a^2-y^2*la+4*z*a*la-z^2*la+4*y*a*la)/a^4;
        dPi_dF.x2222=-1/2*(-2*mu*a^2+4*x*a*la-8*a^2*la-x^2*la-z^2*la+4*z*a*la-2*mu*a^4+2*a^2*la*ln(a^3))/a^4;
        dPi_dF.x3333=-1/2*(-2*mu*a^2+4*x*a*la-8*a^2*la-x^2*la+4*y*a*la-y^2*la-2*mu*a^4+2*a^2*la*ln(a^3))/a^4;
        dPi_dF.x2211=la*(4*a^2-2*a*y+x*y-2*x*a)/a^4;
        dPi_dF.x3311=la*(-2*a*z+x*z-2*x*a+4*a^2)/a^4;
        dPi_dF.x3322=la*(-2*a*z+y*z-2*a*y+4*a^2)/a^4
            dPi_dF.x2121=1/2*(2*x*mu*a^2+2*x*mu*a^4+z^2*x*la-4*x*z*a*la-2*x*y*a*la+8*x*a^2*la-2*x*a^2*la*ln(a^3)-14*a^3*la-2*z^2*a*la+8*z*a^2*la+4*a^3*la*ln(a^3)-4*mu*a^3+2*mu*y*a^2+2*a^4*mu*y+z^2*y*la-4*z*y*a*la+8*y*a^2*la-2*a^2*y*la*ln(a^3))/(x+y)/a^4;
        dPi_dF.x3232=1/2*(2*mu*y*a^2-2*z*y*a*la-4*x*y*a*la+2*a^4*mu*y-2*a^2*y*la*ln(a^3)+8*y*a^2*la+x^2*y*la+8*x*a^2*la-14*a^3*la-4*mu*a^3+4*a^3*la*ln(a^3)-2*x^2*a*la+2*mu*z*a^2-4*x*z*a*la+2*a^4*mu*z-2*z*a^2*la*ln(a^3)+8*z*a^2*la+x^2*z*la)/(y+z)/a^4;
            dPi_dF.x2112=-1/2*(x^2*y*la-2*x^2*a*la+8*x*a^2*la-6*x*y*a*la+y^2*x*la-14*a^3*la-2*z^2*a*la+8*z*a^2*la-2*y^2*a*la+4*a^3*la*ln(a^3)-4*mu*a^3+8*y*a^2*la)/(x+y)/a^4;
        dPi_dF.x3113=-1/2*(-2*x^2*a*la+x^2*z*la+8*x*a^2*la-6*x*z*a*la+z^2*x*la-14*a^3*la-2*z^2*a*la+8*z*a^2*la-2*y^2*a*la+4*a^3*la*ln(a^3)-4*mu*a^3+8*y*a^2*la)/(x+z)/a^4;
        dPi_dF.x3223=-1/2*(-2*y^2*a*la+y^2*z*la+8*y*a^2*la-6*z*y*a*la+z^2*y*la-2*z^2*a*la+8*x*a^2*la-14*a^3*la-4*mu*a^3+4*a^3*la*ln(a^3)-2*x^2*a*la+8*z*a^2*la)/(y+z)/a^4;
            dPi_dF.x3131=1/2*(-2*x*z*a*la+2*x*mu*a^4+y^2*x*la+2*x*mu*a^2-4*x*y*a*la+8*x*a^2*la-2*x*a^2*la*ln(a^3)-4*mu*a^3-14*a^3*la-2*y^2*a*la+8*y*a^2*la+4*a^3*la*ln(a^3)+2*a^4*mu*z+y^2*z*la+2*mu*z*a^2-4*z*y*a*la+8*z*a^2*la-2*z*a^2*la*ln(a^3))/(x+z)/a^4;
    }
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
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
