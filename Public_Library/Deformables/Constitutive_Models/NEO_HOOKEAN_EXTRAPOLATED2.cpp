 //#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Alexey Stomakhin, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/cube.h>
#include <Tools/Math_Tools/pow.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED2.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
NEO_HOOKEAN_EXTRAPOLATED2(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input, const T extra_force_coefficient_input)
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
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
~NEO_HOOKEAN_EXTRAPOLATED2()
{
}
//#####################################################################
// Update Lame Constants
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
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
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    return Energy_Density_Helper(F,id);
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,2>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T x = F.x.x;
    T y = F.x.y;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;
    
    if((dx >= 0) && (dy >= 0))
    {
        T I1=(F*F.Transposed()).Trace(),J=F.Determinant();
        T log_J=log(J);
        return id_mu*((T).5*(I1-TV::m)-log_J)+(T).5*id_lambda*sqr(log_J);
    }
    else if((dx < 0) && (dy >= 0))
    {
        return base.E(a,y) + base.Ex(a,y)*dx + k*dx*dx + id_lambda/a*dx*dx*dy*dy;
    }
    else if((dx >= 0) && (dy < 0))
    {
        return base.E(x,a) + base.Ey(x,a)*dy + k*dy*dy + id_lambda/a*dx*dx*dy*dy;
    }
    else // ((dx < 0) && (dy < 0))
    {
        return base.E(a,a) + base.Ex(a,a)*dx + base.Ey(a,a)*dy + base.Exy(a,a)*dx*dy + k*dx*dx + k*dy*dy + id_lambda/a*dx*dx*dy*dy;
    }
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T x = F.x.x;
    T y = F.x.y;
    T z = F.x.z;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = id_mu;
    T la = id_lambda;   
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;

    
    if((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        T I1=(F*F.Transposed()).Trace(),J=F.Determinant();
        T log_J=log(J);
        return id_mu*((T).5*(I1-TV::m)-log_J)+(T).5*id_lambda*sqr(log_J);
    }
    else if((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
        return (T)0.5*(-mu*cube(a)+mu*a*sqr(y)+mu*a*sqr(z)-mu*a-2*mu*a*log(a*y*z)+la*sqr(log(a*y*z))*a+2*mu*sqr(a)*x-2*mu*x+2*la*log(a*y*z)*x-2*la*log(a*y*z)*a+2*k*a*sqr(x)-4*k*sqr(a)*x+2*k*cube(a))/a;
    }
    else if((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        return (T)0.5*(mu*a*sqr(x)-mu*cube(a)+mu*a*sqr(z)-mu*a-2*mu*a*log(x*a*z)+la*sqr(log(x*a*z))*a+2*y*mu*sqr(a)-2*mu*y+2*la*log(x*a*z)*y-2*la*log(x*a*z)*a+2*sqr(y)*k*a-4*y*k*sqr(a)+2*k*cube(a))/a;
    }
    else if((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        return (T)0.5*(mu*a*sqr(x)+mu*a*sqr(y)-mu*cube(a)-mu*a-2*mu*a*log(x*y*a)+la*sqr(log(x*y*a))*a+2*z*mu*sqr(a)-2*mu*z+2*la*log(x*y*a)*z-2*la*log(x*y*a)*a+2*sqr(z)*k*a-4*z*k*sqr(a)+2*k*cube(a))/a;
    }
    else if((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        return (T)0.5*(mu*sqr(a)-2*x*mu*a-2*x*la*a+2*sqr(y)*k*sqr(a)-2*mu*a*y+sqr(z)*mu*sqr(a)+2*sqr(x)*k*sqr(a)-2*y*la*a+2*a*la*log(sqr(a)*z)*x+2*a*la*log(sqr(a)*z)*y+la*sqr(log(sqr(a)*z))*sqr(a)-2*mu*sqr(a)*log(sqr(a)*z)+2*mu*cube(a)*x-4*la*log(sqr(a)*z)*sqr(a)+2*y*mu*cube(a)+2*la*x*y-4*k*cube(a)*x-4*k*cube(a)*y-2*mu*(a*a*a*a)+2*la*sqr(a)+4*k*(a*a*a*a))/sqr(a);
    }
    else if((dx < 0) && (dy >= 0) && (dz < 0)) // Rzx
    {
        return (T)0.5*(mu*sqr(a)-2*x*mu*a-2*x*la*a+sqr(y)*mu*sqr(a)+2*sqr(z)*k*sqr(a)-2*mu*a*z+2*sqr(x)*k*sqr(a)-2*z*la*a+2*mu*cube(a)*x-4*k*cube(a)*x-2*mu*(a*a*a*a)+2*la*sqr(a)+4*k*(a*a*a*a)+2*z*mu*cube(a)-4*z*k*cube(a)+2*z*la*x+la*sqr(log(sqr(a)*y))*sqr(a)-2*mu*sqr(a)*log(sqr(a)*y)-4*la*log(sqr(a)*y)*sqr(a)+2*a*la*log(sqr(a)*y)*z+2*a*la*log(sqr(a)*y)*x)/sqr(a);
    }
    else if((dx >= 0) && (dy < 0) && (dz < 0)) // Ryz
    {
        return (T)0.5*(mu*sqr(a)+2*sqr(y)*k*sqr(a)-2*mu*a*y+2*sqr(z)*k*sqr(a)-2*mu*a*z+sqr(x)*mu*sqr(a)-2*y*la*a-2*z*la*a+2*y*mu*cube(a)-4*k*cube(a)*y-2*mu*(a*a*a*a)+2*la*sqr(a)+4*k*(a*a*a*a)+2*z*mu*cube(a)+2*z*la*y-4*z*k*cube(a)+2*a*la*log(x*sqr(a))*y+2*a*la*log(x*sqr(a))*z+la*sqr(log(x*sqr(a)))*sqr(a)-2*mu*sqr(a)*log(x*sqr(a))-4*la*log(x*sqr(a))*sqr(a))/sqr(a);
    }
    else // Rxyz
    {
        return (T)0.5*(3*mu*sqr(a)-2*x*mu*a-4*x*la*a+2*sqr(y)*k*sqr(a)-2*mu*a*y+2*sqr(z)*k*sqr(a)-2*mu*a*z+2*sqr(x)*k*sqr(a)-4*y*la*a-4*z*la*a+2*mu*cube(a)*x+2*y*mu*cube(a)+2*la*x*y-4*k*cube(a)*x-4*k*cube(a)*y-3*mu*(a*a*a*a)+6*la*sqr(a)+6*k*(a*a*a*a)+2*z*mu*cube(a)+2*z*la*y-4*z*k*cube(a)+2*z*la*x+2*a*la*log(cube(a))*x+2*a*la*log(cube(a))*y+2*a*la*log(cube(a))*z+la*sqr(log(cube(a)))*sqr(a)-2*mu*sqr(a)*log(cube(a))-6*la*log(cube(a))*sqr(a))/sqr(a);
    }
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    return P_From_Strain_Helper(F,id);
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,2> NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,2>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T x = F.x.x;
    T y = F.x.y;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;
    
    if((dx >= 0) && (dy >= 0))
    {  
        T J=F.Determinant();
        return id_mu*F-(id_mu-id_lambda*log(J))*F.Inverse();
    }
    else if((dx < 0) && (dy >= 0))
    {//[ -(la*(a - s2)^2*(2*a - 2*s1))/a, -(la*(a - s1)^2*(2*a - 2*s2))/a, 0, 0]
        DIAGONAL_MATRIX<T,2> result;
        result.x.x = base.Ex(a,y) + 2*k*dx-id_lambda*(a-y)*(a-y)*(T)2*(a-x)/a;
        result.x.y = base.Ey(a,y) + base.Exy(a,y)*dx - id_lambda*(a-x)*(a-x)*(T)2*(a-y)/a;
        return result;
    }
    else if((dx >= 0) && (dy < 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x.x = base.Ex(x,a) + base.Exy(x,a)*dy -id_lambda*(a-y)*(a-y)*(T)2*(a-x)/a;
        result.x.y = base.Ey(x,a) + 2*k*dy - id_lambda*(a-x)*(a-x)*(T)2*(a-y)/a;
        return result;
    }
    else // ((dx < 0) && (dy < 0))
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x.x = base.Ex(a,a) + base.Exy(a,a)*dy + 2*k*dx -id_lambda*(a-y)*(a-y)*(T)2*(a-x)/a;
        result.x.y = base.Ey(a,a) + base.Exy(a,a)*dx + 2*k*dy - id_lambda*(a-x)*(a-x)*(T)2*(a-y)/a;
        return result;
    }
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T x = F.x.x;
    T y = F.x.y;
    T z = F.x.z;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = id_mu;
    T la = id_lambda;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;
     
    if((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        T J=F.Determinant();
        return id_mu*F-(id_mu-id_lambda*log(J))*F.Inverse();
    }
    else if((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x.x=(mu*sqr(a)-mu+la*log(a*y*z)+2*k*a*x-2*k*sqr(a))/a;
        result.x.y=(mu*a*sqr(y)-mu*a+la*log(a*y*z)*a+la*x-la*a)/a/y;
        result.x.z=(mu*a*sqr(z)-mu*a+la*log(a*y*z)*a+la*x-la*a)/a/z;
        return result;
    }
    else if((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x.x=(mu*a*sqr(x)-mu*a+la*log(x*a*z)*a+la*y-la*a)/x/a;
        result.x.y=(mu*sqr(a)-mu+la*log(x*a*z)+2*k*a*y-2*k*sqr(a))/a;
        result.x.z=(mu*a*sqr(z)-mu*a+la*log(x*a*z)*a+la*y-la*a)/a/z;
        return result;
    }
    else if((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x.x=(mu*a*sqr(x)-mu*a+la*log(x*y*a)*a+la*z-la*a)/x/a;
        result.x.y=(mu*a*sqr(y)-mu*a+la*log(x*y*a)*a+la*z-la*a)/a/y;
        result.x.z=(mu*sqr(a)-mu+la*log(x*y*a)+2*k*a*z-2*k*sqr(a))/a;
        return result;
    }
    else if((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x.x=(mu*cube(a)-mu*a+la*log(sqr(a)*z)*a+la*y-la*a+2*k*sqr(a)*x-2*k*cube(a))/sqr(a);
        result.x.y=(mu*cube(a)-mu*a+la*log(sqr(a)*z)*a+la*x-la*a+2*y*k*sqr(a)-2*k*cube(a))/sqr(a);
        result.x.z=(mu*a*sqr(z)-mu*a+la*log(sqr(a)*z)*a+la*x-2*la*a+la*y)/a/z;
        return result;
    }
    else if((dx < 0) && (dy >= 0) && (dz < 0)) // Rzx
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x.x=(mu*cube(a)-mu*a+la*log(sqr(a)*y)*a+la*z-la*a+2*k*sqr(a)*x-2*k*cube(a))/sqr(a);
        result.x.y=(mu*a*sqr(y)-mu*a+la*log(sqr(a)*y)*a+la*z-2*la*a+la*x)/a/y;
        result.x.z=(mu*cube(a)-mu*a+la*log(sqr(a)*y)*a+la*x-la*a+2*z*k*sqr(a)-2*k*cube(a))/sqr(a);
        return result;
    }
    else if((dx >= 0) && (dy < 0) && (dz < 0)) // Ryz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x.x=(mu*a*sqr(x)-mu*a+la*log(x*sqr(a))*a+la*y-2*la*a+la*z)/x/a;
        result.x.y=(mu*cube(a)-mu*a+la*log(x*sqr(a))*a+la*z-la*a+2*y*k*sqr(a)-2*k*cube(a))/sqr(a);
        result.x.z=(mu*cube(a)-mu*a+la*log(x*sqr(a))*a+la*y-la*a+2*z*k*sqr(a)-2*k*cube(a))/sqr(a);
        return result;
    }
    else // Rxyz
    {
        DIAGONAL_MATRIX<T,3> result;
        result.x.x=(mu*cube(a)-mu*a+la*log(cube(a))*a+la*y-2*la*a+la*z+2*k*sqr(a)*x-2*k*cube(a))/sqr(a);
        result.x.y=(mu*cube(a)-mu*a+la*log(cube(a))*a+la*x-2*la*a+la*z+2*y*k*sqr(a)-2*k*cube(a))/sqr(a);
        result.x.z=(mu*cube(a)-mu*a+la*log(cube(a))*a+la*y-2*la*a+la*x+2*z*k*sqr(a)-2*k*cube(a))/sqr(a);
        return result;
    }
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int id) const
{
    return Isotropic_Stress_Derivative_Helper(F,dP_dF,id);
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T x = F.x.x;
    T y = F.x.y;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;
    
    if((dx >= 0) && (dy >= 0))
    {     
        DIAGONAL_MATRIX<T,2> F_inverse=F.Inverse();
        T mu_minus_lambda_logJ=id_mu+id_lambda*log(F_inverse.Determinant());
        SYMMETRIC_MATRIX<T,2> F_inverse_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(F_inverse.To_Vector());
        dP_dF.x0000=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x00;//alpha+beta+gamma
        dP_dF.x1111=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;
        dP_dF.x1100=id_lambda*F_inverse_outer.x10;//gamma
        dP_dF.x1010=id_mu;//alpha
        dP_dF.x1001=mu_minus_lambda_logJ*F_inverse_outer.x10;//beta
    }
    else if((dx < 0) && (dy >= 0))
    { /* (2*la*(a - s2)^2)/a, (4*la*(a - s1)*(a - s2))/a,                                                     0,                                                     0]
       [ (4*la*(a - s1)*(a - s2))/a,        (2*la*(a - s1)^2)/a,*/
        
        /*   -(2*la*(a - s1)*(a - s2))/(s1 + s2), -(2*la*(a - s1)*(a - s2)*(s1 - a + s2))/(a*(s1 + s2))]*/
        dP_dF.x0000=2*k + (T)2*id_lambda*(a-y)*(a-y)/a;
        dP_dF.x1111=base.Eyy(a,y)+base.Exyy(a,y)*dx + (T)2*id_lambda*(a-x)*(a-x)/a;
        dP_dF.x1100=base.Exy(a,y) + (4*id_lambda*(a - x)*(a -y))/a;
        
        T Ex = base.Ex(a,y)+2*k*dx;
        T Ey = base.Ey(a,y)+base.Exy(a,y)*dx;

        T xpy = x+y; if(fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T xmy = x-y; if(fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;

        dP_dF.x1001=(-Ey*x+Ex*y)/(xpy*xmy) -((T)2*id_lambda*(a - x)*(a - y)*(x - a + y))/(a*(x + y));
        dP_dF.x1010=(-Ey*y+Ex*x)/(xpy*xmy)  -((T)2*id_lambda*(a - x)*(a - y))/(x + y); 
    }
    else if((dx >= 0) && (dy < 0))
    {
        dP_dF.x0000=base.Exx(x,a)+base.Exxy(x,a)*dy+ (T)2*id_lambda*(a-y)*(a-y)/a;
        dP_dF.x1111=2*k+ (T)2*id_lambda*(a-x)*(a-x)/a;
        dP_dF.x1100=base.Exy(x,a)+ (4*id_lambda*(a - x)*(a -y))/a;
        
        T Ex = base.Ex(x,a)+base.Exy(x,a)*dy;
        T Ey = base.Ey(x,a)+2*k*dy;

        T xpy = x+y; if(fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T xmy = x-y; if(fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;

        dP_dF.x1001=(-Ey*x+Ex*y)/(xpy*xmy) -((T)2*id_lambda*(a - x)*(a - y)*(x - a + y))/(a*(x + y));
        dP_dF.x1010=(-Ey*y+Ex*x)/(xpy*xmy)  -((T)2*id_lambda*(a - x)*(a - y))/(x + y); 
    }
    else // ((dx < 0) && (dy < 0))
    {
        dP_dF.x0000=2*k+ (T)2*id_lambda*(a-y)*(a-y)/a;
        dP_dF.x1111=2*k+ (T)2*id_lambda*(a-x)*(a-x)/a;
        dP_dF.x1100=base.Exy(a,a)+ (4*id_lambda*(a - x)*(a -y))/a;

        T xpy = x+y; if(fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
        T cmn = (base.Ex(a,a)-base.Exy(a,a)*a-2*k*a)/xpy;

        dP_dF.x1001=-cmn-base.Exy(a,a) -((T)2*id_lambda*(a - x)*(a - y)*(x - a + y))/(a*(x + y));
        dP_dF.x1010=cmn+2*k-((T)2*id_lambda*(a - x)*(a - y))/(x + y); 
    }

    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    T x = F.x.x;
    T y = F.x.y;
    T z = F.x.z;
    
    T dx = x - extrapolation_cutoff;
    T dy = y - extrapolation_cutoff;
    T dz = z - extrapolation_cutoff;

    T mu = id_mu;
    T la = id_lambda;
    
    T a = extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;

    T xpy = x+y; if(fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
    T xmy = x-y; if(fabs(xmy)<panic_threshold) xmy=xmy<0?-panic_threshold:panic_threshold;
    T xpz = x+z; if(fabs(xpz)<panic_threshold) xpz=xpz<0?-panic_threshold:panic_threshold;
    T xmz = x-z; if(fabs(xmz)<panic_threshold) xmz=xmz<0?-panic_threshold:panic_threshold;
    T ypz = y+z; if(fabs(ypz)<panic_threshold) ypz=ypz<0?-panic_threshold:panic_threshold;
    T ymz = y-z; if(fabs(ymz)<panic_threshold) ymz=ymz<0?-panic_threshold:panic_threshold;

    if((dx >= 0) && (dy >= 0) && (dz >= 0)) // R
    {
        DIAGONAL_MATRIX<T,3> F_inverse=F.Inverse();
        T mu_minus_lambda_logJ=id_mu+id_lambda*log(F_inverse.Determinant());
        SYMMETRIC_MATRIX<T,3> F_inverse_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(F_inverse.To_Vector());
        dP_dF.x0000=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x00;
        dP_dF.x1111=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;
        dP_dF.x2222=id_mu+(id_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
        dP_dF.x1100=id_lambda*F_inverse_outer.x10;
        dP_dF.x2200=id_lambda*F_inverse_outer.x20;
        dP_dF.x2211=id_lambda*F_inverse_outer.x21;
        dP_dF.x1010=id_mu;
        dP_dF.x2020=id_mu;
        dP_dF.x2121=id_mu;
        dP_dF.x1001=mu_minus_lambda_logJ*F_inverse_outer.x10;
        dP_dF.x2002=mu_minus_lambda_logJ*F_inverse_outer.x20;
        dP_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;
    }
    else if((dx < 0) && (dy >= 0) && (dz >= 0)) // Rx
    {
        dP_dF.x0000=2*k;
        dP_dF.x1100=la/y/a;
        dP_dF.x1111=-(-mu*a*sqr(y)-mu*a-2*la*a+la*log(a*y*z)*a+la*x)/sqr(y)/a;
        dP_dF.x2200=la/z/a;
        dP_dF.x2211=la/z/y;
        dP_dF.x2222=-(-mu*a*sqr(z)-mu*a-2*la*a+la*log(a*y*z)*a+la*x)/sqr(z)/a;
        dP_dF.x1001=-(x*mu*a*sqr(y)-x*mu*a+x*la*log(a*y*z)*a+sqr(x)*la-x*la*a-sqr(y)*mu*sqr(a)+mu*sqr(y)-sqr(y)*la*log(a*y*z)-2*sqr(y)*k*a*x+2*sqr(y)*k*sqr(a))/a/y/(xmy*xpy);
        dP_dF.x1010=-(mu*a*sqr(y)-mu*a+la*log(a*y*z)*a+la*x-la*a-mu*sqr(a)*x+mu*x-la*log(a*y*z)*x-2*k*a*sqr(x)+2*k*sqr(a)*x)/a/(xmy*xpy);
        dP_dF.x2002=-(x*mu*a*sqr(z)-x*mu*a+x*la*log(a*y*z)*a+sqr(x)*la-x*la*a-sqr(z)*mu*sqr(a)+mu*sqr(z)-sqr(z)*la*log(a*y*z)-2*sqr(z)*k*a*x+2*sqr(z)*k*sqr(a))/a/z/(xmz*xpz);
        dP_dF.x2020=-(mu*a*sqr(z)-mu*a+la*log(a*y*z)*a+la*x-la*a-mu*sqr(a)*x+mu*x-la*log(a*y*z)*x-2*k*a*sqr(x)+2*k*sqr(a)*x)/a/(xmz*xpz);
        dP_dF.x2112=-(-mu*a+la*log(a*y*z)*a+la*x-la*a)/a/y/z;
        dP_dF.x2121=mu;
    }
    else if((dx >= 0) && (dy < 0) && (dz >= 0)) // Ry
    {
        dP_dF.x0000=-(-mu*a*sqr(x)-mu*a-2*la*a+la*log(x*a*z)*a+la*y)/a/sqr(x);
        dP_dF.x1100=la/x/a;
        dP_dF.x1111=2*k;
        dP_dF.x2200=la/x/z;
        dP_dF.x2211=la/z/a;
        dP_dF.x2222=-(-mu*a*sqr(z)-mu*a-2*la*a+la*log(x*a*z)*a+la*y)/sqr(z)/a;
        dP_dF.x1001=(-sqr(x)*mu*sqr(a)+mu*sqr(x)-sqr(x)*la*log(x*a*z)-2*sqr(x)*k*a*y+2*sqr(x)*k*sqr(a)+y*mu*a*sqr(x)-mu*a*y+y*la*log(x*a*z)*a+sqr(y)*la-y*la*a)/x/a/(xmy*xpy);
        dP_dF.x1010=(-y*mu*sqr(a)+mu*y-la*log(x*a*z)*y-2*sqr(y)*k*a+2*y*k*sqr(a)+mu*a*sqr(x)-mu*a+la*log(x*a*z)*a+la*y-la*a)/a/(xmy*xpy);
        dP_dF.x2002=-(-mu*a+la*log(x*a*z)*a+la*y-la*a)/x/a/z;
        dP_dF.x2020=mu;
        dP_dF.x2112=-(y*mu*a*sqr(z)-mu*a*y+y*la*log(x*a*z)*a+sqr(y)*la-y*la*a-sqr(z)*mu*sqr(a)+mu*sqr(z)-sqr(z)*la*log(x*a*z)-2*sqr(z)*k*a*y+2*sqr(z)*k*sqr(a))/a/z/(ymz*ypz);
        dP_dF.x2121=-(mu*a*sqr(z)-mu*a+la*log(x*a*z)*a+la*y-la*a-y*mu*sqr(a)+mu*y-la*log(x*a*z)*y-2*sqr(y)*k*a+2*y*k*sqr(a))/a/(ymz*ypz);
    }
    else if((dx >= 0) && (dy >= 0) && (dz < 0)) // Rz
    {
        dP_dF.x0000=-(-mu*a*sqr(x)-mu*a-2*la*a+la*log(x*y*a)*a+la*z)/a/sqr(x);
        dP_dF.x1100=la/y/x;
        dP_dF.x1111=-(-mu*a*sqr(y)-mu*a-2*la*a+la*log(x*y*a)*a+la*z)/sqr(y)/a;
        dP_dF.x2200=la/x/a;
        dP_dF.x2211=la/y/a;
        dP_dF.x2222=2*k;
        dP_dF.x1001=-(-mu*a+la*log(x*y*a)*a+la*z-la*a)/x/y/a;
        dP_dF.x1010=mu;
        dP_dF.x2002=(-sqr(x)*mu*sqr(a)+mu*sqr(x)-sqr(x)*la*log(x*y*a)-2*sqr(x)*k*a*z+2*sqr(x)*k*sqr(a)+z*mu*a*sqr(x)-mu*a*z+z*la*log(x*y*a)*a+sqr(z)*la-z*la*a)/x/a/(xmz*xpz);
        dP_dF.x2020=(-z*mu*sqr(a)+mu*z-la*log(x*y*a)*z-2*sqr(z)*k*a+2*z*k*sqr(a)+mu*a*sqr(x)-mu*a+la*log(x*y*a)*a+la*z-la*a)/a/(xmz*xpz);
        dP_dF.x2112=(-sqr(y)*mu*sqr(a)+mu*sqr(y)-sqr(y)*la*log(x*y*a)-2*sqr(y)*k*a*z+2*sqr(y)*k*sqr(a)+z*mu*a*sqr(y)-mu*a*z+z*la*log(x*y*a)*a+sqr(z)*la-z*la*a)/a/y/(ymz*ypz);
        dP_dF.x2121=(-z*mu*sqr(a)+mu*z-la*log(x*y*a)*z-2*sqr(z)*k*a+2*z*k*sqr(a)+mu*a*sqr(y)-mu*a+la*log(x*y*a)*a+la*z-la*a)/a/(ymz*ypz);
        
    }
    else if((dx < 0) && (dy < 0) && (dz >= 0)) // Rxy
    {
        dP_dF.x0000=2*k;
        dP_dF.x1100=la/sqr(a);
        dP_dF.x1111=2*k;
        dP_dF.x2200=la/z/a;
        dP_dF.x2211=la/z/a;
        dP_dF.x2222=-(-mu*a*sqr(z)-mu*a-3*la*a+la*log(sqr(a)*z)*a+la*x+la*y)/sqr(z)/a;
        dP_dF.x1001=-(la*x+mu*cube(a)-mu*a-la*a-2*k*cube(a)+la*log(sqr(a)*z)*a+la*y)/xpy/sqr(a);
        dP_dF.x1010=(2*k*a*x+mu*sqr(a)-mu+la*log(sqr(a)*z)-la-2*k*sqr(a)+2*k*a*y)/xpy/a;
        dP_dF.x2002=-(x*sqr(a)*mu*sqr(z)-mu*sqr(a)*x+x*sqr(a)*la*log(sqr(a)*z)+sqr(x)*la*a-2*x*sqr(a)*la+x*a*la*y-sqr(z)*mu*cube(a)+mu*a*sqr(z)-sqr(z)*la*log(sqr(a)*z)*a-sqr(z)*la*y+sqr(z)*la*a-2*sqr(z)*k*sqr(a)*x+2*sqr(z)*k*cube(a))/sqr(a)/z/(xmz*xpz);
        dP_dF.x2020=-(sqr(z)*mu*sqr(a)-mu*sqr(a)+la*log(sqr(a)*z)*sqr(a)+2*x*la*a-2*la*sqr(a)+y*la*a-mu*cube(a)*x+x*mu*a-a*la*log(sqr(a)*z)*x-la*x*y-2*sqr(x)*k*sqr(a)+2*k*cube(a)*x)/sqr(a)/(xmz*xpz);
        dP_dF.x2112=-(sqr(a)*y*mu*sqr(z)-y*mu*sqr(a)+sqr(a)*y*la*log(sqr(a)*z)+x*a*la*y-2*sqr(a)*y*la+sqr(y)*la*a-sqr(z)*mu*cube(a)+mu*a*sqr(z)-sqr(z)*la*log(sqr(a)*z)*a-sqr(z)*la*x+sqr(z)*la*a-2*sqr(z)*y*k*sqr(a)+2*sqr(z)*k*cube(a))/sqr(a)/z/(ymz*ypz);
        dP_dF.x2121=-(sqr(z)*mu*sqr(a)-mu*sqr(a)+la*log(sqr(a)*z)*sqr(a)+x*la*a-2*la*sqr(a)+2*y*la*a-y*mu*cube(a)+mu*a*y-a*la*log(sqr(a)*z)*y-la*x*y-2*sqr(y)*k*sqr(a)+2*k*cube(a)*y)/sqr(a)/(ymz*ypz);
    }
    else if((dx < 0) && (dy >= 0) && (dz < 0)) // Rzx
    {
        dP_dF.x0000=2*k;
        dP_dF.x1100=la/y/a;
        dP_dF.x1111=-(-mu*a*sqr(y)-mu*a-3*la*a+la*log(sqr(a)*y)*a+la*z+la*x)/sqr(y)/a;
        dP_dF.x2200=la/sqr(a);
        dP_dF.x2211=la/y/a;
        dP_dF.x2222=2*k;
        dP_dF.x1001=-(x*sqr(a)*mu*sqr(y)-mu*sqr(a)*x+x*sqr(a)*la*log(sqr(a)*y)+x*a*la*z-2*x*sqr(a)*la+sqr(x)*la*a-sqr(y)*mu*cube(a)+mu*a*sqr(y)-sqr(y)*la*log(sqr(a)*y)*a-sqr(y)*la*z+sqr(y)*la*a-2*sqr(y)*k*sqr(a)*x+2*sqr(y)*k*cube(a))/sqr(a)/y/(xmy*xpy);
        dP_dF.x1010=-(sqr(y)*mu*sqr(a)-mu*sqr(a)+la*log(sqr(a)*y)*sqr(a)+z*la*a-2*la*sqr(a)+2*x*la*a-mu*cube(a)*x+x*mu*a-a*la*log(sqr(a)*y)*x-z*la*x-2*sqr(x)*k*sqr(a)+2*k*cube(a)*x)/sqr(a)/(xmy*xpy);
        dP_dF.x2002=-(la*x+mu*cube(a)-mu*a-la*a-2*k*cube(a)+la*log(sqr(a)*y)*a+la*z)/xpz/sqr(a);
        dP_dF.x2020=(2*k*a*x+mu*sqr(a)-mu+la*log(sqr(a)*y)-la-2*k*sqr(a)+2*k*a*z)/xpz/a;
        dP_dF.x2112=(-sqr(y)*mu*cube(a)+mu*a*sqr(y)-sqr(y)*la*log(sqr(a)*y)*a-sqr(y)*la*x+sqr(y)*la*a-2*sqr(y)*z*k*sqr(a)+2*sqr(y)*k*cube(a)+sqr(a)*z*mu*sqr(y)-z*mu*sqr(a)+sqr(a)*z*la*log(sqr(a)*y)+sqr(z)*la*a-2*sqr(a)*z*la+x*a*la*z)/sqr(a)/y/(ymz*ypz);
        dP_dF.x2121=(-z*mu*cube(a)+mu*a*z-a*la*log(sqr(a)*y)*z-z*la*x+2*z*la*a-2*sqr(z)*k*sqr(a)+2*z*k*cube(a)+sqr(y)*mu*sqr(a)-mu*sqr(a)+la*log(sqr(a)*y)*sqr(a)-2*la*sqr(a)+x*la*a)/sqr(a)/(ymz*ypz);
    }
    else if((dx >= 0) && (dy < 0) && (dz < 0)) // Ryz
    {
        dP_dF.x0000=-(-mu*a*sqr(x)-mu*a-3*la*a+la*log(x*sqr(a))*a+la*y+la*z)/a/sqr(x);
        dP_dF.x1100=la/x/a;
        dP_dF.x1111=2*k;
        dP_dF.x2200=la/x/a;
        dP_dF.x2211=la/sqr(a);
        dP_dF.x2222=2*k;
        dP_dF.x1001=(-sqr(x)*mu*cube(a)+mu*a*sqr(x)-sqr(x)*la*log(x*sqr(a))*a-sqr(x)*la*z+sqr(x)*la*a-2*sqr(x)*y*k*sqr(a)+2*sqr(x)*k*cube(a)+sqr(a)*y*mu*sqr(x)-y*mu*sqr(a)+sqr(a)*y*la*log(x*sqr(a))+sqr(y)*la*a-2*sqr(a)*y*la+a*y*la*z)/x/sqr(a)/(xmy*xpy);
        dP_dF.x1010=(-y*mu*cube(a)+mu*a*y-a*la*log(x*sqr(a))*y-z*la*y+2*y*la*a-2*sqr(y)*k*sqr(a)+2*k*cube(a)*y+sqr(x)*mu*sqr(a)-mu*sqr(a)+la*log(x*sqr(a))*sqr(a)-2*la*sqr(a)+z*la*a)/sqr(a)/(xmy*xpy);
        dP_dF.x2002=(-sqr(x)*mu*cube(a)+mu*a*sqr(x)-sqr(x)*la*log(x*sqr(a))*a-sqr(x)*la*y+sqr(x)*la*a-2*sqr(x)*z*k*sqr(a)+2*sqr(x)*k*cube(a)+sqr(a)*z*mu*sqr(x)-z*mu*sqr(a)+sqr(a)*z*la*log(x*sqr(a))+a*y*la*z-2*sqr(a)*z*la+sqr(z)*la*a)/x/sqr(a)/(xmz*xpz);
        dP_dF.x2020=(-z*mu*cube(a)+mu*a*z-a*la*log(x*sqr(a))*z-z*la*y+2*z*la*a-2*sqr(z)*k*sqr(a)+2*z*k*cube(a)+sqr(x)*mu*sqr(a)-mu*sqr(a)+la*log(x*sqr(a))*sqr(a)+y*la*a-2*la*sqr(a))/sqr(a)/(xmz*xpz);
        dP_dF.x2112=-(la*y+mu*cube(a)-mu*a-la*a-2*k*cube(a)+la*log(x*sqr(a))*a+la*z)/ypz/sqr(a);
        dP_dF.x2121=(2*k*a*y+mu*sqr(a)-mu+la*log(x*sqr(a))-la-2*k*sqr(a)+2*k*a*z)/ypz/a;
    }
    else // Rxyz
    {
        dP_dF.x0000=2*k;
        dP_dF.x1100=la/sqr(a);
        dP_dF.x1111=2*k;
        dP_dF.x2200=la/sqr(a);
        dP_dF.x2211=la/sqr(a);
        dP_dF.x2222=2*k;
        dP_dF.x1001=-(la*x-mu*a+la*log(cube(a))*a+la*z-2*k*cube(a)-2*la*a+mu*cube(a)+la*y)/xpy/sqr(a);
        dP_dF.x1010=(2*k*sqr(a)*x-mu*a+la*log(cube(a))*a+la*z-2*k*cube(a)-2*la*a+mu*cube(a)+2*y*k*sqr(a))/xpy/sqr(a);
        dP_dF.x2002=-(la*x-mu*a+la*log(cube(a))*a+la*z-2*k*cube(a)-2*la*a+mu*cube(a)+la*y)/xpz/sqr(a);
        dP_dF.x2020=(2*k*sqr(a)*x-mu*a+la*log(cube(a))*a+mu*cube(a)-2*k*cube(a)+la*y-2*la*a+2*z*k*sqr(a))/xpz/sqr(a);
        dP_dF.x2112=-(la*x-mu*a+la*log(cube(a))*a+la*z-2*k*cube(a)-2*la*a+mu*cube(a)+la*y)/ypz/sqr(a);
        dP_dF.x2121=(2*y*k*sqr(a)-mu*a+la*log(cube(a))*a+la*x-2*k*cube(a)-2*la*a+mu*cube(a)+2*z*k*sqr(a))/ypz/sqr(a);
    }
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*id_beta*strain_rate+id_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
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
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED2<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(id_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
namespace PhysBAM{
template class NEO_HOOKEAN_EXTRAPOLATED2<float,2>;
template class NEO_HOOKEAN_EXTRAPOLATED2<float,3>;
template class NEO_HOOKEAN_EXTRAPOLATED2<double,2>;
template class NEO_HOOKEAN_EXTRAPOLATED2<double,3>;
}
