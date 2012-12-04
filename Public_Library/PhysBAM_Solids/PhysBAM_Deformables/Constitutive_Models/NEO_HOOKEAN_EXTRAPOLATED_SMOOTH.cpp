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
NEO_HOOKEAN_EXTRAPOLATED_SMOOTH(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input)
    :youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),
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
        return (T)0.5*(mu*sqr(y)*sqr(a)+mu*sqr(a)-2*mu*sqr(a)*log(y*a)+la*sqr(a)*sqr(log(y*a))-3*sqr(a)*la*log(y*a)-4*mu*x*a+4*a*la*log(y*a)*x+la*sqr(a)-2*x*a*la+sqr(x)*mu*sqr(a)+mu*sqr(x)+sqr(x)*la-la*log(y*a)*sqr(x))/sqr(a);
    }
    else if ((dx >= 0) && (dy < 0))
    {
        return (T)0.5*(sqr(x)*mu*sqr(a)+mu*sqr(a)-2*mu*sqr(a)*log(x*a)+la*sqr(a)*sqr(log(x*a))-3*sqr(a)*la*log(x*a)-4*mu*y*a+4*a*y*la*log(x*a)+la*sqr(a)-2*y*a*la+mu*sqr(y)*sqr(a)+mu*sqr(y)+sqr(y)*la-sqr(y)*la*log(x*a))/sqr(a);
    }
    else // ((dx < 0) && (dy < 0))
    {
        return (T)0.25*(2*mu*sqr(y)*sqr(a)-8*mu*y*cube(a)-8*x*mu*cube(a)+2*sqr(x)*mu*sqr(a)+8*mu*(a*a*a*a)+8*cube(a)*y*la*log(sqr(a))-2*sqr(a)*sqr(y)*la*log(sqr(a))+8*x*cube(a)*la*log(sqr(a))+16*x*y*sqr(a)*la-4*a*sqr(y)*la*x-2*sqr(x)*sqr(a)*la*log(sqr(a))-4*sqr(x)*y*a*la+13*(a*a*a*a)*la+5*sqr(a)*sqr(y)*la+2*(a*a*a*a)*mu*sqr(y)-16*cube(a)*y*la-12*(a*a*a*a)*la*log(sqr(a))-4*(a*a*a*a)*mu*log(sqr(a))+2*(a*a*a*a)*la*sqr(log(sqr(a)))-16*x*la*cube(a)+5*sqr(x)*la*sqr(a)+2*sqr(x)*mu*(a*a*a*a)+sqr(x)*sqr(y)*la)/(a*a*a*a);
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
    T z = F.x33;
    
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
        return (T)0.5*(mu*sqr(y)*sqr(a)+mu*sqr(z)*sqr(a)-2*mu*log(a*y*z)*sqr(a)+la*sqr(log(a*y*z))*sqr(a)-4*a*mu*x+4*a*x*la*log(a*y*z)-3*sqr(a)*la*log(a*y*z)+sqr(x)*mu*sqr(a)+mu*sqr(x)+sqr(x)*la-sqr(x)*la*log(a*y*z)-2*x*a*la+sqr(a)*la)/sqr(a);
    }
    else if ((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        return (T)0.5*(sqr(x)*mu*sqr(a)+mu*sqr(z)*sqr(a)-2*mu*log(x*a*z)*sqr(a)+la*sqr(log(x*a*z))*sqr(a)-4*y*mu*a+4*a*y*la*log(x*a*z)-3*sqr(a)*la*log(x*a*z)+mu*sqr(y)*sqr(a)+mu*sqr(y)+sqr(y)*la-sqr(y)*la*log(x*a*z)-2*y*a*la+sqr(a)*la)/sqr(a);
    }
    else if ((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        return (T)0.5*(sqr(x)*mu*sqr(a)+mu*sqr(y)*sqr(a)-2*mu*log(x*y*a)*sqr(a)+la*sqr(log(x*y*a))*sqr(a)-4*z*mu*a+4*a*z*la*log(x*y*a)-3*sqr(a)*la*log(x*y*a)+mu*sqr(z)*sqr(a)+mu*sqr(z)+sqr(z)*la-sqr(z)*la*log(x*y*a)-2*z*a*la+sqr(a)*la)/sqr(a);
    }
    else if ((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        return (T)0.25*(8*cube(a)*y*la*log(z*sqr(a))-2*sqr(a)*sqr(y)*la*log(z*sqr(a))+8*x*cube(a)*la*log(z*sqr(a))+16*x*y*sqr(a)*la-2*sqr(x)*sqr(a)*la*log(z*sqr(a))+2*mu*sqr(y)*sqr(a)-8*x*mu*cube(a)+2*sqr(x)*mu*sqr(a)+6*mu*(a*a*a*a)-4*sqr(y)*x*a*la+5*sqr(y)*sqr(a)*la+sqr(y)*sqr(x)*la-8*mu*y*cube(a)-4*sqr(x)*y*a*la+5*sqr(x)*sqr(a)*la+13*(a*a*a*a)*la+2*(a*a*a*a)*mu*sqr(z)-4*(a*a*a*a)*mu*log(z*sqr(a))+2*(a*a*a*a)*la*sqr(log(z*sqr(a)))-12*(a*a*a*a)*la*log(z*sqr(a))+2*(a*a*a*a)*mu*sqr(y)-16*cube(a)*y*la-16*x*cube(a)*la+2*sqr(x)*mu*(a*a*a*a))/(a*a*a*a);
    }
    else if ((dx < 0) && (dy >= 0) && (dz < 0)) // Rxz
    {
        return (T)0.25*(2*mu*sqr(z)*sqr(a)-8*x*mu*cube(a)+2*sqr(x)*mu*sqr(a)+6*mu*(a*a*a*a)-4*sqr(z)*x*a*la+5*sqr(z)*sqr(a)*la+sqr(z)*sqr(x)*la+5*sqr(x)*sqr(a)*la-8*mu*z*cube(a)-4*sqr(x)*z*a*la+13*(a*a*a*a)*la+2*(a*a*a*a)*mu*sqr(z)+2*(a*a*a*a)*mu*sqr(y)-16*x*cube(a)*la+2*sqr(x)*mu*(a*a*a*a)-16*z*cube(a)*la+16*z*x*sqr(a)*la-4*(a*a*a*a)*mu*log(y*sqr(a))+2*(a*a*a*a)*la*sqr(log(y*sqr(a)))-12*(a*a*a*a)*la*log(y*sqr(a))+8*cube(a)*z*la*log(y*sqr(a))-2*sqr(a)*sqr(z)*la*log(y*sqr(a))+8*x*cube(a)*la*log(y*sqr(a))-2*sqr(x)*sqr(a)*la*log(y*sqr(a)))/(a*a*a*a);
    }
    else if ((dx >= 0) && (dy < 0) && (dz < 0)) // Rzy
    {
        return (T)0.25*(2*mu*sqr(y)*sqr(a)+2*mu*sqr(z)*sqr(a)+6*mu*(a*a*a*a)+5*sqr(y)*sqr(a)*la+5*sqr(z)*sqr(a)*la-8*mu*y*cube(a)-4*sqr(z)*y*a*la+sqr(z)*sqr(y)*la-8*mu*z*cube(a)-4*sqr(y)*z*a*la+13*(a*a*a*a)*la+2*(a*a*a*a)*mu*sqr(z)+2*(a*a*a*a)*mu*sqr(y)-16*cube(a)*y*la+2*sqr(x)*mu*(a*a*a*a)-16*z*cube(a)*la+16*z*y*sqr(a)*la-4*(a*a*a*a)*mu*log(x*sqr(a))+2*(a*a*a*a)*la*sqr(log(x*sqr(a)))-12*(a*a*a*a)*la*log(x*sqr(a))+8*cube(a)*y*la*log(x*sqr(a))-2*sqr(a)*sqr(y)*la*log(x*sqr(a))+8*z*cube(a)*la*log(x*sqr(a))-2*sqr(z)*sqr(a)*la*log(x*sqr(a)))/(a*a*a*a);
    }
    else // Rxyz
    {
        return (T)0.25*(16*x*y*sqr(a)*la+2*mu*sqr(y)*sqr(a)+2*mu*sqr(z)*sqr(a)-8*x*mu*cube(a)+2*sqr(x)*mu*sqr(a)+12*mu*(a*a*a*a)-4*sqr(y)*x*a*la-4*sqr(z)*x*a*la+8*sqr(y)*sqr(a)*la+sqr(y)*sqr(x)*la+8*sqr(z)*sqr(a)*la+sqr(z)*sqr(x)*la-8*mu*y*cube(a)-4*sqr(x)*y*a*la-4*sqr(z)*y*a*la+8*sqr(x)*sqr(a)*la+sqr(z)*sqr(y)*la-8*mu*z*cube(a)-4*sqr(x)*z*a*la-4*sqr(y)*z*a*la+33*(a*a*a*a)*la+2*(a*a*a*a)*mu*sqr(z)+2*(a*a*a*a)*mu*sqr(y)-28*cube(a)*y*la-28*x*cube(a)*la+2*sqr(x)*mu*(a*a*a*a)-28*z*cube(a)*la+16*z*y*sqr(a)*la+16*z*x*sqr(a)*la-4*(a*a*a*a)*mu*log(cube(a))+2*(a*a*a*a)*la*sqr(log(cube(a)))-18*(a*a*a*a)*la*log(cube(a))+8*cube(a)*y*la*log(cube(a))-2*sqr(a)*sqr(y)*la*log(cube(a))+8*z*cube(a)*la*log(cube(a))-2*sqr(z)*sqr(a)*la*log(cube(a))+8*x*cube(a)*la*log(cube(a))-2*sqr(x)*sqr(a)*la*log(cube(a)))/(a*a*a*a);
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
        result.x22 = (T)0.5*(2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(y*a)-3*la*sqr(a)+4*x*a*la-sqr(x)*la)/y/sqr(a);
        return scale*result;
    }
    else if ((dx >= 0) && (dy < 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = (T)0.5*(2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*a)-3*la*sqr(a)+4*y*a*la-sqr(y)*la)/x/sqr(a);
        result.x22 = (-2*mu*a+2*a*la*log(x*a)-a*la+mu*sqr(a)*y+mu*y+la*y-la*log(x*a)*y)/sqr(a);
        return scale*result;
    }
    else // ((dx < 0) && (dy < 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = (T)0.5*(-4*mu*cube(a)+2*x*mu*sqr(a)+4*cube(a)*la*log(sqr(a))+8*y*la*sqr(a)-2*sqr(y)*a*la-2*x*sqr(a)*la*log(sqr(a))-4*x*y*a*la-8*la*cube(a)+5*x*la*sqr(a)+2*x*mu*(a*a*a*a)+sqr(y)*la*x)/(a*a*a*a);
        result.x22 = (T)0.5*(2*mu*sqr(a)*y-4*mu*cube(a)+4*cube(a)*la*log(sqr(a))-2*sqr(a)*la*log(sqr(a))*y+8*x*la*sqr(a)-4*x*y*a*la-2*sqr(x)*a*la+5*y*la*sqr(a)+2*mu*(a*a*a*a)*y-8*la*cube(a)+sqr(x)*la*y)/(a*a*a*a);
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
    T z = F.x33;
    
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
        result.x11=(-2*mu*a+2*a*la*log(a*y*z)+x*mu*sqr(a)+mu*x+x*la-x*la*log(a*y*z)-a*la)/sqr(a);
        result.x22=(T)0.5*(2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(a*y*z)+4*x*a*la-3*sqr(a)*la-sqr(x)*la)/y/sqr(a);
        result.x33=(T)0.5*(2*mu*sqr(z)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(a*y*z)+4*x*a*la-3*sqr(a)*la-sqr(x)*la)/z/sqr(a);
        return scale*result;
    }
    else if ((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=(T)0.5*(2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*a*z)+4*y*a*la-3*sqr(a)*la-sqr(y)*la)/x/sqr(a);
        result.x22=(-2*mu*a+2*a*la*log(x*a*z)+mu*y*sqr(a)+mu*y+y*la-y*la*log(x*a*z)-a*la)/sqr(a);
        result.x33=(T)0.5*(2*mu*sqr(z)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*a*z)+4*y*a*la-3*sqr(a)*la-sqr(y)*la)/z/sqr(a);
        return scale*result;
    }
    else if ((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=(T)0.5*(2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*y*a)+4*z*a*la-3*sqr(a)*la-sqr(z)*la)/x/sqr(a);
        result.x22=(T)0.5*(2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*y*a)+4*z*a*la-3*sqr(a)*la-sqr(z)*la)/y/sqr(a);
        result.x33=(-2*mu*a+2*a*la*log(x*y*a)+mu*z*sqr(a)+mu*z+z*la-z*la*log(x*y*a)-a*la)/sqr(a);
        return scale*result;
    }
    else if ((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=(T)0.5*(-4*mu*cube(a)+2*x*mu*sqr(a)-2*sqr(y)*a*la+sqr(y)*x*la-4*x*y*a*la+5*x*sqr(a)*la+4*cube(a)*la*log(z*sqr(a))+8*y*sqr(a)*la-2*x*sqr(a)*la*log(z*sqr(a))-8*cube(a)*la+2*x*mu*(a*a*a*a))/(a*a*a*a);
        result.x22=(T)0.5*(2*mu*y*sqr(a)-4*x*y*a*la+5*y*sqr(a)*la+sqr(x)*y*la-4*mu*cube(a)-2*sqr(x)*a*la+4*cube(a)*la*log(z*sqr(a))-2*sqr(a)*y*la*log(z*sqr(a))+8*x*sqr(a)*la+2*(a*a*a*a)*mu*y-8*cube(a)*la)/(a*a*a*a);
        result.x33=(T)0.5/sqr(a)*(4*y*a*la-sqr(y)*la+4*x*a*la-sqr(x)*la+2*mu*sqr(z)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(z*sqr(a))-6*sqr(a)*la)/z;
        return scale*result;
    }
    else if ((dx < 0) && (dy >= 0) && (dz < 0)) // Rxz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=(T)0.5*(-4*mu*cube(a)+2*x*mu*sqr(a)-2*sqr(z)*a*la+sqr(z)*x*la+5*x*sqr(a)*la-4*x*z*a*la-8*cube(a)*la+2*x*mu*(a*a*a*a)+8*z*sqr(a)*la+4*cube(a)*la*log(y*sqr(a))-2*x*sqr(a)*la*log(y*sqr(a)))/(a*a*a*a);
        result.x22=(T)0.5/sqr(a)*(2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(y*sqr(a))-6*sqr(a)*la+4*z*a*la-sqr(z)*la+4*x*a*la-sqr(x)*la)/y;
        result.x33=(T)0.5*(2*mu*z*sqr(a)-4*x*z*a*la+5*z*sqr(a)*la+sqr(x)*z*la-4*mu*cube(a)-2*sqr(x)*a*la+2*(a*a*a*a)*mu*z-8*cube(a)*la+8*x*sqr(a)*la+4*cube(a)*la*log(y*sqr(a))-2*sqr(a)*z*la*log(y*sqr(a)))/(a*a*a*a);
        return scale*result;
    }
    else if ((dx >= 0) && (dy < 0) && (dz < 0)) // Rzy
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=(T)0.5/sqr(a)*(2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*sqr(a))-6*sqr(a)*la+4*y*a*la-sqr(y)*la+4*z*a*la-sqr(z)*la)/x;
        result.x22=(T)0.5*(2*mu*y*sqr(a)+5*y*sqr(a)*la-4*mu*cube(a)-2*sqr(z)*a*la+sqr(z)*y*la-4*z*y*a*la+2*(a*a*a*a)*mu*y-8*cube(a)*la+8*z*sqr(a)*la+4*cube(a)*la*log(x*sqr(a))-2*sqr(a)*y*la*log(x*sqr(a)))/(a*a*a*a);
        result.x33=(T)0.5*(2*mu*z*sqr(a)+5*z*sqr(a)*la-4*z*y*a*la+sqr(y)*z*la-4*mu*cube(a)-2*sqr(y)*a*la+2*(a*a*a*a)*mu*z-8*cube(a)*la+8*y*sqr(a)*la+4*cube(a)*la*log(x*sqr(a))-2*z*sqr(a)*la*log(x*sqr(a)))/(a*a*a*a);
        return scale*result;
    }
    else // Rxyz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x11=(T)0.5*(-4*mu*cube(a)+2*x*mu*sqr(a)-2*sqr(y)*a*la-2*sqr(z)*a*la+sqr(y)*x*la+sqr(z)*x*la-4*x*y*a*la+8*x*sqr(a)*la-4*x*z*a*la+8*y*sqr(a)*la-14*cube(a)*la+2*x*mu*(a*a*a*a)+8*z*sqr(a)*la+4*cube(a)*la*log(cube(a))-2*x*sqr(a)*la*log(cube(a)))/(a*a*a*a);
        result.x22=(T)0.5*(2*mu*y*sqr(a)-4*x*y*a*la+8*y*sqr(a)*la+sqr(x)*y*la-4*mu*cube(a)-2*sqr(x)*a*la-2*sqr(z)*a*la+sqr(z)*y*la-4*z*y*a*la+8*x*sqr(a)*la+2*(a*a*a*a)*mu*y-14*cube(a)*la+8*z*sqr(a)*la+4*cube(a)*la*log(cube(a))-2*sqr(a)*y*la*log(cube(a)))/(a*a*a*a);
        result.x33=(T)0.5*(2*mu*z*sqr(a)-4*x*z*a*la+8*z*sqr(a)*la+sqr(x)*z*la-4*z*y*a*la+sqr(y)*z*la-4*mu*cube(a)-2*sqr(x)*a*la-2*sqr(y)*a*la+2*(a*a*a*a)*mu*z-14*cube(a)*la+8*x*sqr(a)*la+8*y*sqr(a)*la+4*cube(a)*la*log(cube(a))-2*z*sqr(a)*la*log(cube(a)))/(a*a*a*a);
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
        dP_dF.x2222=-(T)0.5*(-2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)-5*la*sqr(a)+2*sqr(a)*la*log(y*a)+4*x*a*la-sqr(x)*la)/sqr(y)/sqr(a);
        dP_dF.x2211=la*(2*a-x)/y/sqr(a);
  
        T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T xmy = x-y; if (fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;

        dP_dF.x2112=-(T)0.5*(-2*x*mu*sqr(a)+2*x*sqr(a)*la*log(y*a)-3*x*la*sqr(a)+4*sqr(x)*a*la-cube(x)*la+4*sqr(y)*mu*a-4*sqr(y)*a*la*log(y*a)+2*sqr(y)*a*la-2*sqr(y)*mu*x-2*sqr(y)*la*x+2*sqr(y)*la*log(y*a)*x)/y/sqr(a)/(sqr(x)-sqr(y));
        dP_dF.x2121=-(T)0.5*(2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(y*a)-3*la*sqr(a)+6*x*a*la-3*sqr(x)*la+4*mu*x*a-4*a*la*log(y*a)*x-2*sqr(x)*mu*sqr(a)-2*mu*sqr(x)+2*la*log(y*a)*sqr(x))/sqr(a)/(sqr(x)-sqr(y));

    }
    else if ((dx >= 0) && (dy < 0))
    {
        dP_dF.x1111=-(T)0.5*(-2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)-5*la*sqr(a)+2*sqr(a)*la*log(x*a)+4*y*a*la-sqr(y)*la)/sqr(x)/sqr(a);
        dP_dF.x2222=-(-mu*sqr(a)-mu-la+la*log(x*a))/sqr(a);
        dP_dF.x2211=la*(2*a-y)/x/sqr(a);
 
        T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T xmy = x-y; if (fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;

        dP_dF.x2112=(T)0.5*(4*sqr(x)*mu*a-4*sqr(x)*a*la*log(x*a)+2*sqr(x)*a*la-2*sqr(x)*mu*y-2*sqr(x)*la*y+2*sqr(x)*la*log(x*a)*y-2*mu*sqr(a)*y+2*y*sqr(a)*la*log(x*a)-3*y*la*sqr(a)+4*sqr(y)*a*la-cube(y)*la)/x/sqr(a)/(sqr(x)-sqr(y));
        dP_dF.x2121=(T)0.5*(4*mu*y*a-4*a*y*la*log(x*a)+6*y*a*la-2*mu*sqr(y)*sqr(a)-2*mu*sqr(y)-3*sqr(y)*la+2*sqr(y)*la*log(x*a)+2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*a)-3*la*sqr(a))/sqr(a)/(sqr(x)-sqr(y));
    }
    else // ((dx < 0) && (dy < 0))
    {
        dP_dF.x1111=-(T)0.5*(-5*la*sqr(a)+2*la*sqr(a)*log(sqr(a))-2*mu*(a*a*a*a)-2*mu*sqr(a)-sqr(y)*la+4*y*a*la)/(a*a*a*a);
        dP_dF.x2222=-(T)0.5*(-2*mu*sqr(a)+2*la*sqr(a)*log(sqr(a))+4*x*a*la-5*la*sqr(a)-2*mu*(a*a*a*a)-sqr(x)*la)/(a*a*a*a);
        dP_dF.x2211=la*(4*sqr(a)-2*y*a-2*x*a+x*y)/(a*a*a*a);

        T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
  
        dP_dF.x2112=-(T)0.5*(sqr(x)*la*y-2*sqr(x)*a*la+8*x*la*sqr(a)-6*x*y*a*la+sqr(y)*la*x+4*cube(a)*la*log(sqr(a))-2*sqr(y)*a*la-4*mu*cube(a)-8*la*cube(a)+8*y*la*sqr(a))/(x+y)/(a*a*a*a);
        dP_dF.x2121=(T)0.5*(2*mu*x*a-2*x*a*la*log(sqr(a))-2*y*la*x+5*x*a*la+2*x*mu*cube(a)-4*mu*sqr(a)-8*la*sqr(a)+4*la*sqr(a)*log(sqr(a))+2*mu*y*a-2*a*y*la*log(sqr(a))+5*y*a*la+2*mu*y*cube(a))/(x+y)/cube(a);
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
    T z = F.x33;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = constant_mu;
    T la = constant_lambda;
    
    T a = extrapolation_cutoff;
    
    T xpy = x+y; if (fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
    T xmy = x-y; if (fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;
    T xpz = x+z; if (fabs(xpz)<panic_threshold) xpz=xpz<0?-panic_threshold:panic_threshold;
    T xmz = x-z; if (fabs(xmz)<panic_threshold) xmz=xmz<0?-panic_threshold:panic_threshold;
    T ypz = y+z; if (fabs(ypz)<panic_threshold) ypz=ypz<0?-panic_threshold:panic_threshold;
    T ymz = y-z; if (fabs(ymz)<panic_threshold) ymz=ymz<0?-panic_threshold:panic_threshold;

    if ((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        DIAGONAL_MATRIX<T,3> F_inverse=F.Inverse();
        T mu_minus_lambda_logJ=constant_mu+constant_lambda*log(F_inverse.Determinant());
        SYMMETRIC_MATRIX<T,3> F_inverse_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(F_inverse.To_Vector());
        dP_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;
        dP_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
        dP_dF.x3333=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x33;
        dP_dF.x2211=constant_lambda*F_inverse_outer.x21;
        dP_dF.x3311=constant_lambda*F_inverse_outer.x31;
        dP_dF.x3322=constant_lambda*F_inverse_outer.x32;
        dP_dF.x2121=constant_mu;
        dP_dF.x3131=constant_mu;
        dP_dF.x3232=constant_mu;
        dP_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;
        dP_dF.x3113=mu_minus_lambda_logJ*F_inverse_outer.x31;
        dP_dF.x3223=mu_minus_lambda_logJ*F_inverse_outer.x32;
    }
    else if ((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
dP_dF.x1111=-(-mu*sqr(a)-mu-la+la*log(a*y*z))/sqr(a);
dP_dF.x2211=la*(2*a-x)/y/sqr(a);
dP_dF.x2222=-(T)0.5*(-2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)-5*sqr(a)*la+2*sqr(a)*la*log(a*y*z)+4*x*a*la-sqr(x)*la)/sqr(y)/sqr(a);
dP_dF.x3311=la*(2*a-x)/z/sqr(a);
dP_dF.x3322=la/y/z;
dP_dF.x3333=-(T)0.5*(-2*mu*sqr(z)*sqr(a)-2*mu*sqr(a)-5*sqr(a)*la+2*sqr(a)*la*log(a*y*z)+4*x*a*la-sqr(x)*la)/sqr(z)/sqr(a);
dP_dF.x2112=-(T)0.5*(-2*x*mu*sqr(a)+2*x*sqr(a)*la*log(a*y*z)+4*sqr(x)*a*la-3*x*sqr(a)*la-cube(x)*la+4*sqr(y)*mu*a-4*sqr(y)*a*la*log(a*y*z)-2*sqr(y)*mu*x-2*sqr(y)*x*la+2*sqr(y)*x*la*log(a*y*z)+2*sqr(y)*a*la)/y/sqr(a)/(xpy*xmy);
dP_dF.x2121=-(T)0.5*(2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(a*y*z)+6*x*a*la-3*sqr(a)*la-3*sqr(x)*la+4*a*mu*x-4*a*x*la*log(a*y*z)-2*sqr(x)*mu*sqr(a)-2*mu*sqr(x)+2*sqr(x)*la*log(a*y*z))/sqr(a)/(xpy*xmy);
dP_dF.x3113=-(T)0.5*(-2*x*mu*sqr(a)+2*x*sqr(a)*la*log(a*y*z)+4*sqr(x)*a*la-3*x*sqr(a)*la-cube(x)*la+4*sqr(z)*mu*a-4*sqr(z)*a*la*log(a*y*z)-2*sqr(z)*mu*x-2*sqr(z)*x*la+2*sqr(z)*x*la*log(a*y*z)+2*sqr(z)*a*la)/z/sqr(a)/(xpz*xmz);
dP_dF.x3131=-(T)0.5*(2*mu*sqr(z)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(a*y*z)+6*x*a*la-3*sqr(a)*la-3*sqr(x)*la+4*a*mu*x-4*a*x*la*log(a*y*z)-2*sqr(x)*mu*sqr(a)-2*mu*sqr(x)+2*sqr(x)*la*log(a*y*z))/sqr(a)/(xpz*xmz);
dP_dF.x3223=-(T)0.5*(-2*mu*sqr(a)+4*x*a*la-3*sqr(a)*la-sqr(x)*la+2*sqr(a)*la*log(a*y*z))/z/sqr(a)/y;
dP_dF.x3232=mu;
    }
    else if ((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
dP_dF.x1111=-(T)0.5*(-2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)-5*sqr(a)*la+2*sqr(a)*la*log(x*a*z)+4*y*a*la-sqr(y)*la)/sqr(x)/sqr(a);
dP_dF.x2211=la*(2*a-y)/x/sqr(a);
dP_dF.x2222=-(-mu*sqr(a)-mu-la+la*log(x*a*z))/sqr(a);
dP_dF.x3311=la/x/z;
dP_dF.x3322=la*(2*a-y)/z/sqr(a);
dP_dF.x3333=-(T)0.5*(-2*mu*sqr(z)*sqr(a)-2*mu*sqr(a)-5*sqr(a)*la+2*sqr(a)*la*log(x*a*z)+4*y*a*la-sqr(y)*la)/sqr(z)/sqr(a);
dP_dF.x2112=(T)0.5*(4*sqr(x)*mu*a-4*sqr(x)*a*la*log(x*a*z)-2*sqr(x)*mu*y-2*sqr(x)*y*la+2*sqr(x)*y*la*log(x*a*z)+2*sqr(x)*a*la-2*mu*y*sqr(a)+2*y*sqr(a)*la*log(x*a*z)+4*sqr(y)*a*la-3*y*sqr(a)*la-cube(y)*la)/x/sqr(a)/(xpy*xmy);
dP_dF.x2121=(T)0.5*(4*y*mu*a-4*a*y*la*log(x*a*z)-2*mu*sqr(y)*sqr(a)-2*mu*sqr(y)-3*sqr(y)*la+2*sqr(y)*la*log(x*a*z)+6*y*a*la+2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*a*z)-3*sqr(a)*la)/sqr(a)/(xpy*xmy);
dP_dF.x3113=-(T)0.5*(-2*mu*sqr(a)+4*y*a*la-3*sqr(a)*la-sqr(y)*la+2*sqr(a)*la*log(x*a*z))/z/sqr(a)/x;
dP_dF.x3131=mu;
dP_dF.x3223=-(T)0.5*(-2*mu*y*sqr(a)+2*y*sqr(a)*la*log(x*a*z)+4*sqr(y)*a*la-3*y*sqr(a)*la-cube(y)*la+4*sqr(z)*mu*a-4*sqr(z)*a*la*log(x*a*z)-2*sqr(z)*mu*y-2*sqr(z)*y*la+2*sqr(z)*y*la*log(x*a*z)+2*sqr(z)*a*la)/z/sqr(a)/(ypz*ymz);
dP_dF.x3232=-(T)0.5*(2*mu*sqr(z)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*a*z)+6*y*a*la-3*sqr(a)*la-3*sqr(y)*la+4*y*mu*a-4*a*y*la*log(x*a*z)-2*mu*sqr(y)*sqr(a)-2*mu*sqr(y)+2*sqr(y)*la*log(x*a*z))/sqr(a)/(ypz*ymz);
    }
    else if ((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
dP_dF.x1111=-(T)0.5*(-2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)-5*sqr(a)*la+2*sqr(a)*la*log(x*y*a)+4*z*a*la-sqr(z)*la)/sqr(x)/sqr(a);
dP_dF.x2211=la/y/x;
dP_dF.x2222=-(T)0.5*(-2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)-5*sqr(a)*la+2*sqr(a)*la*log(x*y*a)+4*z*a*la-sqr(z)*la)/sqr(y)/sqr(a);
dP_dF.x3311=la*(2*a-z)/x/sqr(a);
dP_dF.x3322=la*(2*a-z)/y/sqr(a);
dP_dF.x3333=-(-mu*sqr(a)-mu-la+la*log(x*y*a))/sqr(a);
dP_dF.x2112=-(T)0.5*(-2*mu*sqr(a)+4*z*a*la-3*sqr(a)*la-sqr(z)*la+2*sqr(a)*la*log(x*y*a))/y/sqr(a)/x;
dP_dF.x2121=mu;
dP_dF.x3113=(T)0.5*(4*sqr(x)*mu*a-4*sqr(x)*a*la*log(x*y*a)-2*sqr(x)*mu*z-2*sqr(x)*z*la+2*sqr(x)*z*la*log(x*y*a)+2*sqr(x)*a*la-2*mu*z*sqr(a)+2*z*sqr(a)*la*log(x*y*a)+4*sqr(z)*a*la-3*z*sqr(a)*la-cube(z)*la)/x/sqr(a)/(xpz*xmz);
dP_dF.x3131=(T)0.5*(4*z*mu*a-4*a*z*la*log(x*y*a)-2*mu*sqr(z)*sqr(a)-2*mu*sqr(z)-3*sqr(z)*la+2*sqr(z)*la*log(x*y*a)+6*z*a*la+2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*y*a)-3*sqr(a)*la)/sqr(a)/(xpz*xmz);
dP_dF.x3223=(T)0.5*(4*sqr(y)*mu*a-4*sqr(y)*a*la*log(x*y*a)-2*sqr(y)*mu*z-2*sqr(y)*z*la+2*sqr(y)*z*la*log(x*y*a)+2*sqr(y)*a*la-2*mu*z*sqr(a)+2*z*sqr(a)*la*log(x*y*a)+4*sqr(z)*a*la-3*z*sqr(a)*la-cube(z)*la)/y/sqr(a)/(ypz*ymz);
dP_dF.x3232=(T)0.5*(4*z*mu*a-4*a*z*la*log(x*y*a)-2*mu*sqr(z)*sqr(a)-2*mu*sqr(z)-3*sqr(z)*la+2*sqr(z)*la*log(x*y*a)+6*z*a*la+2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)+2*sqr(a)*la*log(x*y*a)-3*sqr(a)*la)/sqr(a)/(ypz*ymz);
    }
    else if ((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
dP_dF.x1111=-(T)0.5*(-5*sqr(a)*la+2*sqr(a)*la*log(z*sqr(a))-2*mu*(a*a*a*a)-2*mu*sqr(a)-sqr(y)*la+4*y*a*la)/(a*a*a*a);
dP_dF.x2211=la*(-2*a*y+x*y-2*x*a+4*sqr(a))/(a*a*a*a);
dP_dF.x2222=-(T)0.5*(-2*mu*sqr(a)+4*x*a*la-5*sqr(a)*la-sqr(x)*la+2*sqr(a)*la*log(z*sqr(a))-2*mu*(a*a*a*a))/(a*a*a*a);
dP_dF.x3311=la*(2*a-x)/z/sqr(a);
dP_dF.x3322=la*(2*a-y)/z/sqr(a);
dP_dF.x3333=-(T)0.5/sqr(a)*(4*y*a*la-sqr(y)*la+4*x*a*la-sqr(x)*la-2*mu*sqr(z)*sqr(a)-2*mu*sqr(a)-8*sqr(a)*la+2*sqr(a)*la*log(z*sqr(a)))/sqr(z);
dP_dF.x2112=-(T)0.5*(-2*sqr(x)*a*la+sqr(x)*y*la+8*x*sqr(a)*la-6*x*y*a*la+sqr(y)*x*la-4*mu*cube(a)-8*cube(a)*la-2*sqr(y)*a*la+4*cube(a)*la*log(z*sqr(a))+8*y*sqr(a)*la)/xpy/(a*a*a*a);
dP_dF.x2121=(T)0.5*(5*x*a*la-2*y*x*la-2*x*a*la*log(z*sqr(a))+2*a*mu*x+2*x*mu*cube(a)+4*sqr(a)*la*log(z*sqr(a))-4*mu*sqr(a)-8*sqr(a)*la+5*y*a*la-2*a*y*la*log(z*sqr(a))+2*y*mu*a+2*mu*y*cube(a))/xpy/cube(a);
dP_dF.x3113=-(T)0.5*(4*x*cube(a)*y*la-x*sqr(a)*sqr(y)*la+4*sqr(x)*cube(a)*la-cube(x)*sqr(a)*la-2*x*mu*(a*a*a*a)+2*x*(a*a*a*a)*la*log(z*sqr(a))-6*x*(a*a*a*a)*la+4*sqr(z)*mu*cube(a)-2*x*mu*sqr(z)*sqr(a)+2*sqr(z)*sqr(y)*a*la-sqr(z)*sqr(y)*x*la+4*sqr(z)*x*y*a*la-5*sqr(z)*x*sqr(a)*la-4*sqr(z)*cube(a)*la*log(z*sqr(a))-8*sqr(z)*y*sqr(a)*la+2*sqr(z)*x*sqr(a)*la*log(z*sqr(a))+8*sqr(z)*cube(a)*la)/z/(a*a*a*a)/(xpz*xmz);
dP_dF.x3131=-(T)0.5*(4*cube(a)*y*la-sqr(y)*sqr(a)*la+12*x*cube(a)*la-6*sqr(x)*sqr(a)*la+2*(a*a*a*a)*mu*sqr(z)-2*mu*(a*a*a*a)+2*(a*a*a*a)*la*log(z*sqr(a))-6*(a*a*a*a)*la+4*x*mu*cube(a)-2*sqr(x)*mu*sqr(a)+2*sqr(y)*x*a*la-sqr(y)*sqr(x)*la+4*sqr(x)*y*a*la-4*x*cube(a)*la*log(z*sqr(a))-8*x*y*sqr(a)*la+2*sqr(x)*sqr(a)*la*log(z*sqr(a))-2*sqr(x)*mu*(a*a*a*a))/(a*a*a*a)/(xpz*xmz);
dP_dF.x3223=-(T)0.5*(4*sqr(y)*cube(a)*la-cube(y)*sqr(a)*la+4*x*cube(a)*y*la-y*sqr(a)*sqr(x)*la-2*(a*a*a*a)*mu*y+2*y*(a*a*a*a)*la*log(z*sqr(a))-6*y*(a*a*a*a)*la-2*y*mu*sqr(z)*sqr(a)+4*sqr(z)*x*y*a*la-5*sqr(z)*y*sqr(a)*la-sqr(z)*sqr(x)*y*la+4*sqr(z)*mu*cube(a)+2*sqr(z)*sqr(x)*a*la-4*sqr(z)*cube(a)*la*log(z*sqr(a))+2*sqr(z)*sqr(a)*y*la*log(z*sqr(a))-8*sqr(z)*x*sqr(a)*la+8*sqr(z)*cube(a)*la)/z/(a*a*a*a)/(ypz*ymz);
dP_dF.x3232=-(T)0.5*(12*cube(a)*y*la-6*sqr(y)*sqr(a)*la+4*x*cube(a)*la-sqr(x)*sqr(a)*la+2*(a*a*a*a)*mu*sqr(z)-2*mu*(a*a*a*a)+2*(a*a*a*a)*la*log(z*sqr(a))-6*(a*a*a*a)*la-2*mu*sqr(y)*sqr(a)+4*sqr(y)*x*a*la-sqr(y)*sqr(x)*la+4*mu*y*cube(a)+2*sqr(x)*y*a*la-4*cube(a)*y*la*log(z*sqr(a))+2*sqr(a)*sqr(y)*la*log(z*sqr(a))-8*x*y*sqr(a)*la-2*(a*a*a*a)*mu*sqr(y))/(a*a*a*a)/(ypz*ymz);
    }
    else if ((dx < 0) && (dy >= 0) && (dz < 0)) // Rxz
    {
dP_dF.x1111=-(T)0.5*(-5*sqr(a)*la+2*sqr(a)*la*log(y*sqr(a))-2*mu*(a*a*a*a)-2*mu*sqr(a)-sqr(z)*la+4*z*a*la)/(a*a*a*a);
dP_dF.x2211=la*(2*a-x)/y/sqr(a);
dP_dF.x2222=-(T)0.5/sqr(a)*(-2*mu*sqr(y)*sqr(a)-2*mu*sqr(a)-8*sqr(a)*la+2*sqr(a)*la*log(y*sqr(a))+4*z*a*la-sqr(z)*la+4*x*a*la-sqr(x)*la)/sqr(y);
dP_dF.x3311=la*(-2*a*z+x*z-2*x*a+4*sqr(a))/(a*a*a*a);
dP_dF.x3322=la*(2*a-z)/y/sqr(a);
dP_dF.x3333=-(T)0.5*(-2*mu*sqr(a)+4*x*a*la-5*sqr(a)*la-sqr(x)*la-2*mu*(a*a*a*a)+2*sqr(a)*la*log(y*sqr(a)))/(a*a*a*a);
dP_dF.x2112=-(T)0.5*(-2*x*mu*(a*a*a*a)+2*x*(a*a*a*a)*la*log(y*sqr(a))-6*x*(a*a*a*a)*la+4*x*cube(a)*z*la-sqr(z)*x*sqr(a)*la+4*sqr(x)*cube(a)*la-cube(x)*sqr(a)*la+4*cube(a)*mu*sqr(y)-2*x*mu*sqr(y)*sqr(a)+2*sqr(z)*sqr(y)*a*la-sqr(z)*sqr(y)*x*la-5*x*sqr(a)*sqr(y)*la+4*sqr(y)*x*z*a*la+8*sqr(y)*cube(a)*la-8*sqr(y)*z*sqr(a)*la-4*sqr(y)*cube(a)*la*log(y*sqr(a))+2*sqr(y)*x*sqr(a)*la*log(y*sqr(a)))/y/(a*a*a*a)/(xpy*xmy);
dP_dF.x2121=-(T)0.5*(2*(a*a*a*a)*mu*sqr(y)-2*mu*(a*a*a*a)+2*(a*a*a*a)*la*log(y*sqr(a))-6*(a*a*a*a)*la+4*z*cube(a)*la-sqr(z)*sqr(a)*la+12*x*cube(a)*la-6*sqr(x)*sqr(a)*la+4*x*mu*cube(a)-2*sqr(x)*mu*sqr(a)+2*sqr(z)*x*a*la-sqr(z)*sqr(x)*la+4*sqr(x)*z*a*la-2*sqr(x)*mu*(a*a*a*a)-8*z*x*sqr(a)*la-4*x*cube(a)*la*log(y*sqr(a))+2*sqr(x)*sqr(a)*la*log(y*sqr(a)))/(a*a*a*a)/(xpy*xmy);
dP_dF.x3113=-(T)0.5*(-2*sqr(x)*a*la+sqr(x)*z*la+8*x*sqr(a)*la-6*x*z*a*la+sqr(z)*x*la-4*mu*cube(a)+4*cube(a)*la*log(y*sqr(a))-2*sqr(z)*a*la-8*cube(a)*la+8*z*sqr(a)*la)/xpz/(a*a*a*a);
dP_dF.x3131=(T)0.5*(5*x*a*la-2*z*x*la+2*x*mu*cube(a)+2*a*mu*x-2*x*a*la*log(y*sqr(a))-8*sqr(a)*la-4*mu*sqr(a)+4*sqr(a)*la*log(y*sqr(a))+5*z*a*la+2*mu*z*cube(a)+2*z*mu*a-2*a*z*la*log(y*sqr(a)))/xpz/cube(a);
dP_dF.x3223=(T)0.5*(-2*sqr(y)*mu*z*sqr(a)+4*sqr(y)*x*z*a*la-5*sqr(y)*z*sqr(a)*la-sqr(x)*sqr(y)*z*la+4*cube(a)*mu*sqr(y)+2*sqr(x)*sqr(y)*a*la+8*sqr(y)*cube(a)*la-8*x*sqr(a)*sqr(y)*la-4*sqr(y)*cube(a)*la*log(y*sqr(a))+2*sqr(y)*sqr(a)*z*la*log(y*sqr(a))-2*(a*a*a*a)*mu*z+2*z*(a*a*a*a)*la*log(y*sqr(a))-6*z*(a*a*a*a)*la+4*sqr(z)*cube(a)*la-cube(z)*sqr(a)*la+4*x*cube(a)*z*la-sqr(x)*z*sqr(a)*la)/y/(a*a*a*a)/(ypz*ymz);
dP_dF.x3232=(T)0.5*(-2*mu*sqr(z)*sqr(a)+4*sqr(z)*x*a*la-6*sqr(z)*sqr(a)*la-sqr(z)*sqr(x)*la+4*mu*z*cube(a)+2*sqr(x)*z*a*la-2*(a*a*a*a)*mu*sqr(z)+12*z*cube(a)*la-8*z*x*sqr(a)*la-4*cube(a)*z*la*log(y*sqr(a))+2*sqr(a)*sqr(z)*la*log(y*sqr(a))+2*(a*a*a*a)*mu*sqr(y)-2*mu*(a*a*a*a)+2*(a*a*a*a)*la*log(y*sqr(a))-6*(a*a*a*a)*la+4*x*cube(a)*la-sqr(x)*sqr(a)*la)/(a*a*a*a)/(ypz*ymz);
    }
    else if ((dx >= 0) && (dy < 0) && (dz < 0)) // Rzy
    {
dP_dF.x1111=-(T)0.5/sqr(a)*(-2*sqr(x)*mu*sqr(a)-2*mu*sqr(a)-8*sqr(a)*la+2*sqr(a)*la*log(x*sqr(a))+4*y*a*la-sqr(y)*la+4*z*a*la-sqr(z)*la)/sqr(x);
dP_dF.x2211=la*(2*a-y)/x/sqr(a);
dP_dF.x2222=-(T)0.5*(-2*mu*sqr(a)-5*sqr(a)*la-sqr(z)*la+4*z*a*la-2*mu*(a*a*a*a)+2*sqr(a)*la*log(x*sqr(a)))/(a*a*a*a);
dP_dF.x3311=la*(2*a-z)/x/sqr(a);
dP_dF.x3322=la*(-2*a*z+y*z-2*a*y+4*sqr(a))/(a*a*a*a);
dP_dF.x3333=-(T)0.5*(-5*sqr(a)*la+2*sqr(a)*la*log(x*sqr(a))-2*mu*(a*a*a*a)-2*mu*sqr(a)-sqr(y)*la+4*y*a*la)/(a*a*a*a);
dP_dF.x2112=(T)0.5*(-2*sqr(x)*mu*y*sqr(a)-5*y*sqr(a)*sqr(x)*la+4*sqr(x)*mu*cube(a)+2*sqr(z)*sqr(x)*a*la-sqr(z)*sqr(x)*y*la+4*sqr(x)*z*y*a*la+8*sqr(x)*cube(a)*la-8*sqr(x)*z*sqr(a)*la-4*sqr(x)*cube(a)*la*log(x*sqr(a))+2*sqr(x)*sqr(a)*y*la*log(x*sqr(a))-2*(a*a*a*a)*mu*y+2*y*(a*a*a*a)*la*log(x*sqr(a))-6*y*(a*a*a*a)*la+4*sqr(y)*cube(a)*la-cube(y)*sqr(a)*la+4*y*cube(a)*z*la-sqr(z)*y*sqr(a)*la)/x/(a*a*a*a)/(xpy*xmy);
dP_dF.x2121=(T)0.5*(-2*mu*sqr(y)*sqr(a)-6*sqr(y)*sqr(a)*la+4*mu*y*cube(a)+2*sqr(z)*y*a*la-sqr(z)*sqr(y)*la+4*sqr(y)*z*a*la-2*(a*a*a*a)*mu*sqr(y)+12*cube(a)*y*la-8*z*y*sqr(a)*la-4*cube(a)*y*la*log(x*sqr(a))+2*sqr(a)*sqr(y)*la*log(x*sqr(a))+2*sqr(x)*mu*(a*a*a*a)-2*mu*(a*a*a*a)+2*(a*a*a*a)*la*log(x*sqr(a))-6*(a*a*a*a)*la+4*z*cube(a)*la-sqr(z)*sqr(a)*la)/(a*a*a*a)/(xpy*xmy);
dP_dF.x3113=(T)0.5*(-2*sqr(x)*mu*z*sqr(a)-5*sqr(x)*z*sqr(a)*la+4*sqr(x)*z*y*a*la-sqr(x)*sqr(y)*z*la+4*sqr(x)*mu*cube(a)+2*sqr(x)*sqr(y)*a*la+8*sqr(x)*cube(a)*la-8*y*sqr(a)*sqr(x)*la-4*sqr(x)*cube(a)*la*log(x*sqr(a))+2*sqr(x)*z*sqr(a)*la*log(x*sqr(a))-2*(a*a*a*a)*mu*z+2*z*(a*a*a*a)*la*log(x*sqr(a))-6*z*(a*a*a*a)*la+4*y*cube(a)*z*la-sqr(y)*z*sqr(a)*la+4*sqr(z)*cube(a)*la-cube(z)*sqr(a)*la)/x/(a*a*a*a)/(xpz*xmz);
dP_dF.x3131=(T)0.5*(-2*mu*sqr(z)*sqr(a)-6*sqr(z)*sqr(a)*la+4*sqr(z)*y*a*la-sqr(z)*sqr(y)*la+4*mu*z*cube(a)+2*sqr(y)*z*a*la-2*(a*a*a*a)*mu*sqr(z)+12*z*cube(a)*la-8*z*y*sqr(a)*la-4*z*cube(a)*la*log(x*sqr(a))+2*sqr(z)*sqr(a)*la*log(x*sqr(a))+2*sqr(x)*mu*(a*a*a*a)-2*mu*(a*a*a*a)+2*(a*a*a*a)*la*log(x*sqr(a))-6*(a*a*a*a)*la+4*cube(a)*y*la-sqr(y)*sqr(a)*la)/(a*a*a*a)/(xpz*xmz);
dP_dF.x3223=-(T)0.5*(-2*sqr(y)*a*la+sqr(y)*z*la+8*y*sqr(a)*la-6*z*y*a*la+sqr(z)*y*la-4*mu*cube(a)+4*cube(a)*la*log(x*sqr(a))-2*sqr(z)*a*la-8*cube(a)*la+8*z*sqr(a)*la)/ypz/(a*a*a*a);
dP_dF.x3232=(T)0.5*(-2*z*y*la+2*mu*y*cube(a)+2*y*mu*a+5*y*a*la-2*a*y*la*log(x*sqr(a))-8*sqr(a)*la+4*sqr(a)*la*log(x*sqr(a))-4*mu*sqr(a)+2*mu*z*cube(a)+2*z*mu*a+5*z*a*la-2*z*a*la*log(x*sqr(a)))/ypz/cube(a);
    }
    else // Rxyz
    {
dP_dF.x1111=-(T)0.5*(-8*sqr(a)*la+2*sqr(a)*la*log(cube(a))-2*mu*(a*a*a*a)-2*mu*sqr(a)-sqr(y)*la+4*z*a*la-sqr(z)*la+4*y*a*la)/(a*a*a*a);
dP_dF.x2211=la*(-2*a*y+x*y-2*x*a+4*sqr(a))/(a*a*a*a);
dP_dF.x2222=-(T)0.5*(-2*mu*sqr(a)+4*x*a*la-8*sqr(a)*la-sqr(x)*la-sqr(z)*la+4*z*a*la-2*mu*(a*a*a*a)+2*sqr(a)*la*log(cube(a)))/(a*a*a*a);
dP_dF.x3311=la*(-2*a*z+x*z-2*x*a+4*sqr(a))/(a*a*a*a);
dP_dF.x3322=la*(-2*a*z+y*z-2*a*y+4*sqr(a))/(a*a*a*a);
dP_dF.x3333=-(T)0.5*(-2*mu*sqr(a)+4*x*a*la-8*sqr(a)*la-sqr(x)*la+4*y*a*la-sqr(y)*la-2*mu*(a*a*a*a)+2*sqr(a)*la*log(cube(a)))/(a*a*a*a);
dP_dF.x2112=-(T)0.5*(-2*sqr(x)*a*la+sqr(x)*y*la+8*x*sqr(a)*la-6*x*y*a*la+sqr(y)*x*la-2*sqr(y)*a*la-14*cube(a)*la-2*sqr(z)*a*la+8*z*sqr(a)*la+4*cube(a)*la*log(cube(a))-4*mu*cube(a)+8*y*sqr(a)*la)/xpy/(a*a*a*a);
dP_dF.x2121=(T)0.5*(2*x*mu*sqr(a)+2*x*mu*(a*a*a*a)+sqr(z)*x*la-4*x*z*a*la-2*x*y*a*la+8*x*sqr(a)*la-2*x*sqr(a)*la*log(cube(a))-14*cube(a)*la-2*sqr(z)*a*la+8*z*sqr(a)*la+4*cube(a)*la*log(cube(a))-4*mu*cube(a)+2*mu*y*sqr(a)+2*(a*a*a*a)*mu*y+sqr(z)*y*la-4*z*y*a*la+8*y*sqr(a)*la-2*sqr(a)*y*la*log(cube(a)))/xpy/(a*a*a*a);
dP_dF.x3113=-(T)0.5*(-2*sqr(x)*a*la+sqr(x)*z*la+8*x*sqr(a)*la-6*x*z*a*la+sqr(z)*x*la-2*sqr(y)*a*la-14*cube(a)*la-2*sqr(z)*a*la+8*z*sqr(a)*la+4*cube(a)*la*log(cube(a))-4*mu*cube(a)+8*y*sqr(a)*la)/xpz/(a*a*a*a);
dP_dF.x3131=(T)0.5*(-2*x*z*a*la+sqr(y)*x*la+2*x*mu*sqr(a)-4*x*y*a*la+2*x*mu*(a*a*a*a)+8*x*sqr(a)*la-2*x*sqr(a)*la*log(cube(a))-2*sqr(y)*a*la-14*cube(a)*la+4*cube(a)*la*log(cube(a))-4*mu*cube(a)+8*y*sqr(a)*la+sqr(y)*z*la+2*mu*z*sqr(a)-4*z*y*a*la+2*(a*a*a*a)*mu*z+8*z*sqr(a)*la-2*z*sqr(a)*la*log(cube(a)))/xpz/(a*a*a*a);
dP_dF.x3223=-(T)0.5*(-2*sqr(y)*a*la+sqr(y)*z*la+8*y*sqr(a)*la-6*z*y*a*la+sqr(z)*y*la-2*sqr(x)*a*la+8*x*sqr(a)*la-14*cube(a)*la-4*mu*cube(a)+4*cube(a)*la*log(cube(a))-2*sqr(z)*a*la+8*z*sqr(a)*la)/ypz/(a*a*a*a);
dP_dF.x3232=(T)0.5*(-2*z*y*a*la+2*mu*y*sqr(a)+8*y*sqr(a)*la+2*(a*a*a*a)*mu*y-2*sqr(a)*y*la*log(cube(a))+sqr(x)*y*la-4*x*y*a*la-2*sqr(x)*a*la+8*x*sqr(a)*la-14*cube(a)*la-4*mu*cube(a)+4*cube(a)*la*log(cube(a))+2*mu*z*sqr(a)+8*z*sqr(a)*la+2*(a*a*a*a)*mu*z-2*z*sqr(a)*la*log(cube(a))+sqr(x)*z*la-4*x*z*a*la)/ypz/(a*a*a*a);
    }
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
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
namespace PhysBAM{
template class NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<float,2>;
template class NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<float,3>;
template class NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<double,2>;
template class NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<double,3>;
}
