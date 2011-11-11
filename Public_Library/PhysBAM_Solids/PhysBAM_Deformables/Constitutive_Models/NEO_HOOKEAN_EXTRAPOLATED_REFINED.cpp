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
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_REFINED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
NEO_HOOKEAN_EXTRAPOLATED_REFINED(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input,const T corner_cutoff_input, const T extra_force_coefficient_input):
    youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),
    extrapolation_cutoff(extrapolation_cutoff_input),corner_cutoff(corner_cutoff_input),extra_force_coefficient(extra_force_coefficient_input),
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
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
~NEO_HOOKEAN_EXTRAPOLATED_REFINED()
{
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    return Energy_Density_Helper(F,simplex);
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,2>& F,const int simplex) const
{
    T x = F.x11;
    T y = F.x22;
    
    T a = extrapolation_cutoff;
    T r = corner_cutoff;
    T k = extra_force_coefficient;

    if ((x>=a) && (y>=a) && (x+y>=r+2*a)) // main region
    {    
        return base.E(x,y);
    }
    else if ((y>=a+r) && (x<a)) // top
    {
        T dx = x-a;
        return base.E(a,y) + base.Ex(a,y)*dx + sqr(dx)*(base.Exx(a,a+r)-base.Eyy(a,a+r))/4 + k*sqr(dx);
    }
    else if ((x>=a+r) && (y<a)) // right
    {                                 
        T dy = y-a;            
        return base.E(x,a) + base.Ey(x,a)*dy -sqr(dy)*(base.Exx(a+r,a)-base.Eyy(a+r,a))/4 + k*sqr(dy);
    }            
    else if ((y<=r+x) && (x<=r+y) && (x+y<r+2*a)) // center
    {
        T dx = (x+y-r-2*a)/2;
        T dy = (x+y-r-2*a)/2;
        T x0 = x-dx;
        T y0 = y-dy;
        return base.E(x0,y0) + base.Ex(x0,y0)*dx + base.Ey(x0,y0)*dy + base.Exy(a,a+r)*dx*dy + 2*k*dx*dy;
    }          
    else if ((y>x+r) && (y<a+r)) // corner top
    {       
        T dx = x-a;
        T dy = y-a-r;
        return base.E(a,a+r) + base.Ex(a,a+r)*dx + base.Ey(a,a+r)*dy + base.Exy(a,a+r)*dx*dy + (sqr(dx)-sqr(dy))*(base.Exx(a,a+r)-base.Eyy(a,a+r))/4 + k*(sqr(dx)+sqr(dy));
    }       
    else if ((x>y+r) && (x<a+r)) // corner right
    {   
        T dx = x-a-r;
        T dy = y-a;
        return base.E(a+r,a) + base.Ex(a+r,a)*dx + base.Ey(a+r,a)*dy + base.Exy(a+r,a)*dx*dy + (sqr(dx)-sqr(dy))*(base.Exx(a+r,a)-base.Eyy(a+r,a))/4 + k*(sqr(dx)+sqr(dy));
    }
    else
    {
        PHYSBAM_FATAL_ERROR();
        return 0;
    }
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return 0;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    return P_From_Strain_Helper(F,scale,simplex);
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,2> NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,2>& F,const T scale,const int simplex) const
{
    T x = F.x11;
    T y = F.x22;
    
    T a = extrapolation_cutoff;
    T r = corner_cutoff;
    T k = extra_force_coefficient;

    T mu = constant_mu;
    T la = constant_lambda;

    T la_log_J = la*log(x*y);

    if ((x>=a) && (y>=a) && (x+y>=r+2*a)) // main region
    {    
        T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda,J=F.Determinant();
        return scale_mu*F-(scale_mu-scale_lambda*log(J))*F.Inverse();
    }
    else if ((y>=a+r) && (x<a)) // top
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = (mu * x * x - mu + la_log_J + 2 * k * x * x - 2 * k * x * a) / x;
        result.x22 = (mu * y * y - mu + la_log_J) / y;
        return scale*result;
    }
    else if ((x>=a+r) && (y<a)) // right
    {                                 
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = (mu * x * x - mu + la_log_J) / x;
        result.x22 = (mu * y * y - mu + la_log_J + 2 * k * y * y - 2 * k * y * a) / y;
        return scale*result;
    }            
    else if ((y<=r+x) && (x<=r+y) && (x+y<r+2*a)) // center
    {
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = (mu * x * x - mu + la_log_J + k * x * x + k * x * y - k * x * r - 2 * k * x * a) / x;
        result.x22 = (mu * y * y - mu + la_log_J + k * x * y + k * y * y - k * y * r - 2 * k * y * a) / y;
        return scale*result;
    }          
    else if ((y>x+r) && (y<a+r)) // corner top
    {       
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = (mu * x * x - mu + la_log_J + 2 * k * x * x - 2 * k * x * a) / x;
        result.x22 = (mu * y * y - mu + la_log_J + 2 * k * y * y - 2 * k * y * a - 2 * k * y * r) / y;
        return scale*result;
    }       
    else if ((x>y+r) && (x<a+r)) // corner right
    {   
        DIAGONAL_MATRIX<T,2> result;
        result.x11 = (mu * x * x - mu + la_log_J + 2 * k * x * x - 2 * k * x * a - 2 * k * x * r) / x;
        result.x22 = (mu * y * y - mu + la_log_J + 2 * k * y * y - 2 * k * y * a) / y;
        return scale*result;
    }
    else
    {
        PHYSBAM_FATAL_ERROR();
        return DIAGONAL_MATRIX<T,2>();
    }
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return DIAGONAL_MATRIX<T,3>();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int triangle) const
{
    return Isotropic_Stress_Derivative_Helper(F,dP_dF,triangle);
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    T x = F.x11;
    T y = F.x22;
    
    T a = extrapolation_cutoff;
    T r = corner_cutoff;
    T k = extra_force_coefficient;

    T mu = constant_mu;
    T la = constant_lambda;

    //T la_log_J = la*log(x*y);

    if ((x>=a) && (y>=a) && (x+y>=r+2*a)) // main region
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
    else if ((y>=a+r) && (x<a)) // top
    {
        dP_dF.x1111 = -(-mu * x * x - mu - la + la * log(x * y) - 0.2e1 * k * x * x) * pow(x, -0.2e1);
        dP_dF.x2222 = -(-mu * y * y - mu - la + la * log(x * y)) * pow(y, -0.2e1);
        dP_dF.x2211 = la / y / x;
        dP_dF.x2121 = -(mu * y * y - mu * x * x - 0.2e1 * k * x * x + 0.2e1 * k * x * a) / (x * x - y * y);
        dP_dF.x2112 = -(-mu * x * x + x * x * la * log(x * y) + mu * y * y - y * y * la * log(x * y) - 0.2e1 * y * y * k * x * x + 0.2e1 * y * y * k * x * a) / x / y / (x * x - y * y);
    }
    else if ((x>=a+r) && (y<a)) // right
    {                                 
        dP_dF.x1111 = -(-mu * x * x - mu - la + la * log(x * y)) * pow(x, -0.2e1);
        dP_dF.x2222 = -(-mu * y * y - mu - la + la * log(x * y) - 0.2e1 * k * y * y) * pow(y, -0.2e1);
        dP_dF.x2211 = la / y / x;
        dP_dF.x2121 = (-mu * y * y - 0.2e1 * k * y * y + 0.2e1 * k * y * a + mu * x * x) / (x * x - y * y);
        dP_dF.x2112 = -(-mu * x * x + x * x * la * log(x * y) + 0.2e1 * y * y * k * x * x - 0.2e1 * x * x * k * y * a + mu * y * y - y * y * la * log(x * y)) / x / y / (x * x - y * y);
    }            
    else if ((y<=r+x) && (x<=r+y) && (x+y<r+2*a)) // center
    {
        dP_dF.x1111 = -(-mu * x * x - mu - la + la * log(x * y) - k * x * x) * pow(x, -0.2e1);
        dP_dF.x2222 = -(-mu * y * y - mu - la + la * log(x * y) - k * y * y) * pow(y, -0.2e1);
        dP_dF.x2211 = (la + k * x * y) / x / y;
        dP_dF.x2121 = -(-x * mu - k * x + k * r + 0.2e1 * k * a - y * mu - k * y) / (x + y);
        dP_dF.x2112 = -(x * x * k * y - x * mu - x * k * y * r - 0.2e1 * x * k * y * a + x * la * log(x * y) + y * y * k * x - y * mu + y * la * log(x * y)) / (x + y) / x / y;
    }          
    else if ((y>x+r) && (y<a+r)) // corner top
    {       
        dP_dF.x1111 = -(-mu * x * x - mu - la + la * log(x * y) - 0.2e1 * k * x * x) * pow(x, -0.2e1);
        dP_dF.x2222 = -(-mu * y * y - mu - la + la * log(x * y) - 0.2e1 * k * y * y) * pow(y, -0.2e1);
        dP_dF.x2211 = la / y / x;
        dP_dF.x2121 = -(mu * y * y + 0.2e1 * k * y * y - 0.2e1 * k * y * a - 0.2e1 * k * y * r - mu * x * x - 0.2e1 * k * x * x + 0.2e1 * k * x * a) / (x * x - y * y);
        dP_dF.x2112 = -(-mu * x * x + x * x * la * log(x * y) - 0.2e1 * x * x * k * y * a - 0.2e1 * x * x * k * y * r + mu * y * y - y * y * la * log(x * y) + 0.2e1 * y * y * k * x * a) / x / y / (x * x - y * y);
    }       
    else if ((x>y+r) && (x<a+r)) // corner right
    {   
        dP_dF.x1111 = -(-mu * x * x - mu - la + la * log(x * y) - 0.2e1 * k * x * x) * pow(x, -0.2e1);
        dP_dF.x2222 = -(-mu * y * y - mu - la + la * log(x * y) - 0.2e1 * k * y * y) * pow(y, -0.2e1);
        dP_dF.x2211 = la / y / x;
        dP_dF.x2121 = -(mu * y * y + 0.2e1 * k * y * y - 0.2e1 * k * y * a - mu * x * x - 0.2e1 * k * x * x + 0.2e1 * k * x * a + 0.2e1 * k * x * r) / (x * x - y * y);
        dP_dF.x2112 = -(-mu * x * x + x * x * la * log(x * y) - 0.2e1 * x * x * k * y * a + mu * y * y - y * y * la * log(x * y) + 0.2e1 * y * y * k * x * a + 0.2e1 * y * y * k * x * r) / x / y / (x * x - y * y);
    }
    else
    {
        PHYSBAM_FATAL_ERROR();
    }

    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int triangle) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
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
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
template class NEO_HOOKEAN_EXTRAPOLATED_REFINED<float,2>;
template class NEO_HOOKEAN_EXTRAPOLATED_REFINED<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NEO_HOOKEAN_EXTRAPOLATED_REFINED<double,2>;
template class NEO_HOOKEAN_EXTRAPOLATED_REFINED<double,3>;
#endif
