//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Russell Howes, Joseph Teran.
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
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input, const T extra_force_coefficient_input):
youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),
extrapolation_cutoff(extrapolation_cutoff_input), extra_force_coefficient(extra_force_coefficient_input),
panic_threshold((T)1e-6)
{
    assert(poissons_ratio>-1&&poissons_ratio<.5);
    constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    constant_mu=youngs_modulus/(2*(1+poissons_ratio));
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;
    base.Initialize(constant_mu,constant_lambda);
    //s1o=-1; s2o=-1; s3o=-1;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
~NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA()
{
}
//#####################################################################
// Update Lame Constants
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::Update_Lame_Constants(const T youngs_modulus_input, const T poissons_ratio_input,const T Rayleigh_coefficient_input)
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
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    T I1=(F*F.Transposed()).Trace(),J=F.Determinant();
    T c = extrapolation_cutoff;    
    
    if (J>=c){
        T log_J=log(J);
        return constant_mu*((T).5*(I1-TV::m)-log_J)+(T).5*constant_lambda*sqr(log_J);
    }
    return Energy_Density_Helper(F,simplex);
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,2>& F,const int simplex) const
{
    
    T c = extrapolation_cutoff;    
    T s1 = F.x11;
    T s2 = F.x22;
    T mu = constant_mu;
    T la = constant_lambda;
    
    
    T t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16;//,t17,t18,t19,t0;
    T t0;
    
    t2 = s1-s2;
    t3 = s2/2;
    t4 = c*4;
    t5 = t2*t2;
    t6 = t4+t5;
    t7 = sqrt(t6);
    t8 = t7/2;
    t10 = s1/2;
    t9 = -t10+t3+t8;
    t11 = 1/t9;
    t12 = t10-t3+t8;
    t13 = t12*t9;
    t14 = log(t13);
    t15 = 1/t12;
    t16 = t10+t3-t8;
    t0 = -mu*t14+t16*(mu*t12-mu*t15+la*t14*t15)+t16*(-mu*t11+mu*t9+la*t11*t14)+la*(t14*t14)/2+mu*(t12*t12+t9*t9-2)/2;
    
    
    return t0;
    
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int simplex) const
{
    T a,a1,a2,a3,a11,a12,a13,a22,a23,a33;    
    T s1 = F.x11;
    T s2 = F.x22;
    T s3 = F.x33;
    /////////////Hopefully I can get the member function to work at some point
    int maxiter = 100; int iter;
    T tol = (T)1e-8;
    //if (d==2){PHYSBAM_FATAL_ERROR();}
    //if (pow(s1-s1o,2)+pow(s1-s1o,2)+pow(s1-s1o,2)<tol) {return;} //Don't need to recalculate
    T c = extrapolation_cutoff;    
    T c13 = pow(c,(T)1/3);
    T a0,f,fa;
    
    if (s1<0) {a0=c13-s1;}          //Pick initial guess
    else{ if (s2<0) {a0=c13-s2;}
    else{ if (s3<0) {a0=c13-s3;}
        else a0=c13;}}
    
    iter=0; f=(T)2e-8; a=a0;
    
    while(f>tol){                   //Newton iteration
        f=(s1+a)*(s2+a)*(s3+a)-c;
        fa=(s1+a)*(s2+a)+(s1+a)*(s3+a)+(s3+a)*(s2+a);
        if (abs(fa)<tol){PHYSBAM_FATAL_ERROR();}
        a=a-f/fa;
        iter++; if (iter>maxiter){PHYSBAM_FATAL_ERROR();}
    }
    
    //Now we get a1,a2,a3,a11,a12,a13,a22,a23,a33
    
    T den = (3*a*a+2*a*(s1+s2+s3)+(s1*s2+s1*s3+s2*s3));
    
    a1 = -(a+s2)*(a+s3)/den;
    a2 = -(a+s1)*(a+s3)/den;
    a3 = -(a+s1)*(a+s2)/den;
    a11 = -2*a1*(s3+2*a+a1*s3+3*a1*a+s2+a1*s2+a1*s1)/den;
    a22 = -2*a2*(s3+2*a+a2*s3+3*a2*a+s1+a2*s1+a2*s2)/den;
    a33 = -2*a3*(s1+2*a+a3*s1+3*a3*a+s2+a3*s2+a3*s3)/den;
    a12 = -(s3+a+a1*s3+2*a1*a+a1*s1+2*a1*a2*s3+6*a1*a2*a+2*a1*a2*s2+2*a1*a2*s1+a2*s3+2*a2*a+a2*s2)/den;
    a13 = -(s2+a+a1*s2+2*a1*a+a1*s1+2*a1*a3*s2+6*a1*a3*a+2*a1*a3*s3+2*a1*a3*s1+a3*s2+2*a3*a+a3*s3)/den;
    a23 = -(s1+a+a3*s1+2*a3*a+a3*s3+2*a3*a2*s1+6*a3*a2*a+2*a3*a2*s2+2*a3*a2*s3+a2*s1+2*a2*a+a2*s2)/den;
    
    T mu = constant_mu;
    T la = constant_lambda;
    //Calculate_A(s1,s2,s3);
    T k = extra_force_coefficient*youngs_modulus; 
    
    if(simplex==-2){ //This is a signal to print out a and its derivatives
        printf("a  : %10.5e \n",a);
        printf("a1 : %10.5e \n",a1);
        printf("a2 : %10.5e \n",a2);
        printf("a3 : %10.5e \n",a3);
        printf("a11: %10.5e \n",a11);
        printf("a12: %10.5e \n",a12);
        printf("a13: %10.5e \n",a13);
        printf("a22: %10.5e \n",a22);
        printf("a23: %10.5e \n",a23);
        printf("a33: %10.5e \n",a33);
        printf("s1 : %10.5e \n",s1);
        printf("s2 : %10.5e \n",s2);
        printf("s3 : %10.5e \n",s3);    
        printf("mu : %10.5e \n",mu);
        printf("la : %10.5e \n",la);
        printf("c  : %10.5e \n",c);    
        printf("k  : %10.5e \n",k);
    }    
    /////////////Hopefully I can get the member function to work at some point   
    

    
    T t0 = -mu*log((a+s1)*(a+s2)*(a+s3))-a*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+(a*a)*k+la*pow(log((a+s1)*(a+s2)*(a+s3)),2)/2+mu*(pow(a+s1,2)+pow(a+s2,2)+pow(a+s3,2)-2)/2;
    
    return t0;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda,J=F.Determinant();
    if(J>=extrapolation_cutoff) return scale_mu*F-(scale_mu-scale_lambda*log(J))*F.Inverse();
    return P_From_Strain_Helper(F,scale,simplex);
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,2> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,2>& F,const T scale,const int simplex) const
{
    T s1 = F.x11;
    T s2 = F.x22;
    
    T c = extrapolation_cutoff;
    
    T mu = constant_mu;
    T la = constant_lambda;
    DIAGONAL_MATRIX<T,2> result;
    T t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67;//,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87;
    
    t18 = s1-s2;
    t19 = c*4;
    t20 = t18*t18;
    t21 = t19+t20;
    t22 = s1*2;
    t23 = s2*2;
    t24 = t22-t23;
    t25 = 1/sqrt(t21);
    t26 = t24*t25/4;
    t27 = t26-(T)0.5;
    t28 = s1/2;
    t29 = s2/2;
    t30 = sqrt(t21);
    t31 = t30/2;
    t32 = -t28+t29+t31;
    t33 = 1/(t32*t32);
    t34 = t28-t29+t31;
    t35 = t26+(T)0.5;
    t36 = t27*t34;
    t37 = t32*t35;
    t38 = t36+t37;
    t39 = 1/(t34*t34);
    t40 = t32*t34;
    t41 = log(t40);
    t42 = t24*t25/2;
    t43 = 1/t32;
    t44 = 1/t34;
    t45 = t28+t29-t31;
    t46 = mu*t27;
    t47 = mu*t27*t33;
    t48 = la*t33*t38*t44;
    t49 = t46+t47+t48-la*t27*t33*t41;
    t50 = mu*t35;
    t51 = mu*t35*t39;
    t52 = la*t38*t39*t43;
    t53 = t50+t51+t52-la*t35*t39*t41;
    t54 = t45*t53;
    t55 = t42-1;
    t56 = t32*t55;
    t57 = t42+1;
    t58 = t34*t57;
    t59 = t56+t58;
    t60 = mu*t59/2;
    t61 = mu*t32;
    t62 = la*t41*t43;
    t63 = t61+t62-mu*t43;
    t64 = mu*t34;
    t65 = la*t41*t44;
    t66 = t64+t65-mu*t44;
    t67 = la*t38*t41*t43*t44;
    result.x11 = t54+t60+t67-t27*t63-t27*t66+t49*(t28+t29-t30/2)-mu*t38*t43*t44;
    result.x22 = -t54-t60-t67-t45*t49+t35*t63+t35*t66+mu*t38*t43*t44;
    
    
    return scale*result;
    
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const
{
    T s1 = F.x11;
    T s2 = F.x22;
    T s3 = F.x33;
    
    /////////////Hopefully I can get the member function to work at some point
    T a,a1,a2,a3,a11,a12,a13,a22,a23,a33;
    int maxiter = 100; int iter;
    T tol = (T)1e-8;
    T c = extrapolation_cutoff;    
    T c13 = pow(c,(T)1/3);
    T a0,f,fa;
    
    if (s1<0) {a0=c13-s1;}          //Pick initial guess
    else{ if (s2<0) {a0=c13-s2;}
    else{ if (s3<0) {a0=c13-s3;}
        else a0=c13;}}
    
    iter=0; f=(T)2e-8; a=a0;
    
    while(f>tol){                   //Newton iteration
        f=(s1+a)*(s2+a)*(s3+a)-c;
        fa=(s1+a)*(s2+a)+(s1+a)*(s3+a)+(s3+a)*(s2+a);
        if (abs(fa)<tol){PHYSBAM_FATAL_ERROR();}
        a=a-f/fa;
        iter++; if (iter>maxiter){PHYSBAM_FATAL_ERROR();}
    }
    
    //Now we get a1,a2,a3,a11,a12,a13,a22,a23,a33
    
    T den = (3*a*a+2*a*(s1+s2+s3)+(s1*s2+s1*s3+s2*s3));
    
    a1 = -(a+s2)*(a+s3)/den;
    a2 = -(a+s1)*(a+s3)/den;
    a3 = -(a+s1)*(a+s2)/den;
    a11 = -2*a1*(s3+2*a+a1*s3+3*a1*a+s2+a1*s2+a1*s1)/den;
    a22 = -2*a2*(s3+2*a+a2*s3+3*a2*a+s1+a2*s1+a2*s2)/den;
    a33 = -2*a3*(s1+2*a+a3*s1+3*a3*a+s2+a3*s2+a3*s3)/den;
    a12 = -(s3+a+a1*s3+2*a1*a+a1*s1+2*a1*a2*s3+6*a1*a2*a+2*a1*a2*s2+2*a1*a2*s1+a2*s3+2*a2*a+a2*s2)/den;
    a13 = -(s2+a+a1*s2+2*a1*a+a1*s1+2*a1*a3*s2+6*a1*a3*a+2*a1*a3*s3+2*a1*a3*s1+a3*s2+2*a3*a+a3*s3)/den;
    a23 = -(s1+a+a3*s1+2*a3*a+a3*s3+2*a3*a2*s1+6*a3*a2*a+2*a3*a2*s2+2*a3*a2*s3+a2*s1+2*a2*a+a2*s2)/den;
    
    /////////////Hopefully I can get the member function to work at some point
    
    T mu = constant_mu;
    T la = constant_lambda;
    DIAGONAL_MATRIX<T,3> result;
    //Calculate_A(s1,s2,s3);
    T k = extra_force_coefficient*youngs_modulus;    
    result.x11 = -a*(mu*(a1+1)+mu*1/pow(a+s1,2)*(a1+1)-la*log((a+s1)*(a+s2)*(a+s3))*1/pow(a+s1,2)*(a1+1)+(la*1/pow(a+s1,2)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1)))/((a+s2)*(a+s3)))-a1*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a1*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a1*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+mu*((a1+1)*(a*2+s1*2)+a1*(a*2+s2*2)+a1*(a*2+s3*2))/2-a*(a1*mu+a1*mu*1/pow(a+s2,2)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1/pow(a+s2,2)+(la*1/pow(a+s2,2)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1)))/((a+s1)*(a+s3)))-a*(a1*mu+a1*mu*1/pow(a+s3,2)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1/pow(a+s3,2)+(la*1/pow(a+s3,2)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1)))/((a+s1)*(a+s2)))+a*a1*k*2-(mu*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1)))/((a+s1)*(a+s2)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1)))/((a+s1)*(a+s2)*(a+s3));
    result.x22 = -a*(mu*(a2+1)+mu*1/pow(a+s2,2)*(a2+1)-la*log((a+s1)*(a+s2)*(a+s3))*1/pow(a+s2,2)*(a2+1)+(la*1/pow(a+s2,2)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1)))/((a+s1)*(a+s3)))-a2*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a2*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a2*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+mu*((a2+1)*(a*2+s2*2)+a2*(a*2+s1*2)+a2*(a*2+s3*2))/2-a*(a2*mu+a2*mu*1/pow(a+s1,2)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1/pow(a+s1,2)+(la*1/pow(a+s1,2)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1)))/((a+s2)*(a+s3)))-a*(a2*mu+a2*mu*1/pow(a+s3,2)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1/pow(a+s3,2)+(la*1/pow(a+s3,2)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1)))/((a+s1)*(a+s2)))+a*a2*k*2-(mu*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1)))/((a+s1)*(a+s2)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1)))/((a+s1)*(a+s2)*(a+s3));
    result.x33 = -a*(mu*(a3+1)+mu*1/pow(a+s3,2)*(a3+1)-la*log((a+s1)*(a+s2)*(a+s3))*1/pow(a+s3,2)*(a3+1)+(la*1/pow(a+s3,2)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1)))/((a+s1)*(a+s2)))-a3*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a3*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a3*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+mu*((a3+1)*(a*2+s3*2)+a3*(a*2+s1*2)+a3*(a*2+s2*2))/2-a*(a3*mu+a3*mu*1/pow(a+s1,2)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1/pow(a+s1,2)+(la*1/pow(a+s1,2)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1)))/((a+s2)*(a+s3)))-a*(a3*mu+a3*mu*1/pow(a+s2,2)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1/pow(a+s2,2)+(la*1/pow(a+s2,2)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1)))/((a+s1)*(a+s3)))+a*a3*k*2-(mu*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1)))/((a+s1)*(a+s2)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1)))/((a+s1)*(a+s2)*(a+s3));
    
    return scale*result;
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int triangle) const
{
    
    Isotropic_Stress_Derivative_Helper(F,dP_dF,triangle);    
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    
    T c = extrapolation_cutoff;
    T mu = constant_mu;
    T la = constant_lambda;
    T J=F.Determinant();
    
    if(J>=c){
        DIAGONAL_MATRIX<T,2> F_inverse=F.Inverse();
        T mu_minus_lambda_logJ=constant_mu-constant_lambda*log(J);
        SYMMETRIC_MATRIX<T,2> F_inverse_outer=SYMMETRIC_MATRIX<T,2>::Outer_Product(F_inverse.To_Vector());
        dP_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;//alpha+beta+gamma
        dP_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
        dP_dF.x2211=constant_lambda*F_inverse_outer.x21;//gamma
        dP_dF.x2121=constant_mu;//alpha
        dP_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;//beta
        return;
    }
    
    T s1 = F.x11;
    T s2 = F.x22;
    
    T t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115;
    
    
    t69 = s1-s2;
    t70 = c*4;
    t71 = t69*t69;
    t72 = t70+t71;
    t73 = pow(t72,(T)1.5);
    t74 = s1*s1;
    t75 = s2*s2;
    t76 = t74*t74;
    t77 = t75*t75;
    t78 = sqrt(t72);
    t79 = s2/2;
    t80 = t78/2;
    t81 = s1/2;
    t82 = -t79+t80+t81;
    t83 = t79+t80-t81;
    t84 = t82*t83;
    t85 = log(t84);
    t86 = t74*t75*12;
    t87 = s1*t73*2;
    t88 = c*t74*16;
    t89 = c*t75*16;
    t90 = c*c;
    t91 = t90*32;
    t92 = t76*2;
    t93 = t77*2;
    t94 = c*mu*32;
    t95 = mu*t76*2;
    t96 = mu*t77*2;
    t97 = mu*t74*t75*12;
    t98 = la*s1*s2*t85*16;
    t99 = la*s1*t78*t85*4;
    t100 = la*s2*t78*t85*4;
    t101 = mu*s1*t75*t78;
    t102 = mu*s2*t74*t78;
    t104 = s1*s2*2;
    t103 = -t104+t70+t74+t75;
    t105 = pow(t103,(T)1.5);
    t106 = log(c);
    t107 = 1/c;
    t108 = sqrt(t103);
    t109 = s1+s2;if(fabs(t109)<panic_threshold) t109=t109<0?-panic_threshold:panic_threshold;
    t110 = 1/t109;
    t111 = 1/sqrt(t103);
    t112 = mu*t90*2;
    t113 = mu*s1*t108;
    t114 = mu*s2*t108;
    t115 = c*la*t106*2;
    dP_dF.x1111 = (t100+t101+t102+t94+t95+t96+t97+t98+t99+mu*t74*4+mu*t75*12-la*t74*log(t82*(-s1/2+t79+t80))*4-c*la*t85*32+c*mu*t74*12+c*mu*t75*4-mu*s1*s2*16-la*t75*t85*12+mu*s1*t73*3-mu*s1*t78*4-mu*s2*t73-mu*s2*t78*4-c*mu*s1*s2*16-mu*s1*s2*t74*8-mu*s1*s2*t75*8-mu*s1*t74*t78-mu*s2*t75*t78)/(t86+t87+t88+t89+t91+t92+t93-s2*t73*2-c*s1*s2*32-s1*s2*t74*8-s1*s2*t75*8);
    dP_dF.x2222 = (t100+t101+t102+t94+t95+t96+t97+t98+t99+mu*t74*12+mu*t75*4-c*la*t85*32+c*mu*t74*4+c*mu*t75*12-mu*s1*s2*16-la*t74*t85*12-la*t75*t85*4-mu*s1*t73-mu*s1*t78*4+mu*s2*t73*3-mu*s2*t78*4-c*mu*s1*s2*16-mu*s1*s2*t74*8-mu*s1*s2*t75*8-mu*s1*t74*t78-mu*s2*t75*t78)/(t86-t87+t88+t89+t91+t92+t93+s2*t73*2-c*s1*s2*32-s1*s2*t74*8-s1*s2*t75*8);
    dP_dF.x2211 = -1/pow(t103,(T)1.5)*t107*(mu*t105-c*mu*s1*2-c*mu*s2*2-la*t105*t106+mu*s1*t90*2+mu*s2*t90*2+c*la*s1*t106*2+c*la*s2*t106*2);
    dP_dF.x2121 = t107*t110*t111*(t112+t113+t114+t115-c*mu*2-mu*t74-mu*t75+c*mu*t74+c*mu*t75+la*t106*t74+la*t106*t75-la*s1*t106*t108-la*s2*t106*t108);
    dP_dF.x2112 = -t107*t110*t111*(t112-t113-t114+t115-c*mu*2+mu*s1*s2*2-c*mu*s1*s2*2-la*s1*s2*t106*2+la*s1*t106*t108+la*s2*t106*t108);
    
    
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int triangle) const
{
    
    T c = extrapolation_cutoff;
    T mu = constant_mu;
    T la = constant_lambda;
    T J=F.Determinant();
    T k = extra_force_coefficient*youngs_modulus;
    
    if(J>=c){
        DIAGONAL_MATRIX<T,3> F_inverse=F.Inverse();
        T mu_minus_lambda_logJ=constant_mu+constant_lambda*log(F_inverse.Determinant());
        SYMMETRIC_MATRIX<T,3> F_inverse_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(F_inverse.To_Vector());
        dP_dF.x1111=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x11;
        dP_dF.x2222=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x22;
        dP_dF.x3333=constant_mu+(constant_lambda+mu_minus_lambda_logJ)*F_inverse_outer.x33;
        dP_dF.x2211=constant_lambda*F_inverse_outer.x21;
        dP_dF.x3311=constant_lambda*F_inverse_outer.x31;
        dP_dF.x3322=constant_lambda*F_inverse_outer.x32;
        dP_dF.x2121=constant_mu;dP_dF.x3131=constant_mu;dP_dF.x3232=constant_mu;
        dP_dF.x2112=mu_minus_lambda_logJ*F_inverse_outer.x21;
        dP_dF.x3113=mu_minus_lambda_logJ*F_inverse_outer.x31;
        dP_dF.x3223=mu_minus_lambda_logJ*F_inverse_outer.x32;
        return;
    }
    
    T s1 = F.x11;
    T s2 = F.x22;
    T s3 = F.x33;
    
    //Now we need to check and make sure these aren't too close to each other or their negatives
    
    
    /////////////Hopefully I can get the member function to work at some point
    T a,a1,a2,a3,a11,a12,a13,a22,a23,a33;
    int maxiter = 100; int iter;
    T tol = (T)1e-8;
    T c13 = pow(c,(T)1/3);
    T a0,f,fa;
    
    if (s1<0) {a0=c13-s1;}          //Pick initial guess
    else{ if (s2<0) {a0=c13-s2;}
    else{ if (s3<0) {a0=c13-s3;}
        else a0=c13;}}
    
    iter=0; f=(T)2e-8; a=a0;
    
    while(f>tol){                   //Newton iteration
        f=(s1+a)*(s2+a)*(s3+a)-c;
        fa=(s1+a)*(s2+a)+(s1+a)*(s3+a)+(s3+a)*(s2+a);
        if (abs(fa)<tol){PHYSBAM_FATAL_ERROR();}
        a=a-f/fa;
        iter++; if (iter>maxiter){PHYSBAM_FATAL_ERROR();}
    }
    
    //This part is super jimmy-rigged
    

    
    T s1ms3 = s1-s3; if(fabs(s1ms3)<panic_threshold) {s1=s1+1*panic_threshold; s3=s3+3*panic_threshold; s1ms3=s1-s3;}
    T s2ms3 = s2-s3; if(fabs(s2ms3)<panic_threshold) {s2=s2+2*panic_threshold; s3=s3+3*panic_threshold; s2ms3=s2-s3;}
    T s1ps2 = s1+s2; if(fabs(s1ps2)<panic_threshold) {s1=s1+1*panic_threshold; s2=s2+2*panic_threshold; s1ps2=s1+s2;}
    T s1ps3 = s1+s3; if(fabs(s1ps3)<panic_threshold) {s1=s1+1*panic_threshold; s3=s3+3*panic_threshold; s1ps3=s1+s3;}
    T s2ps3 = s2+s3; if(fabs(s2ps3)<panic_threshold) {s2=s2+2*panic_threshold; s3=s3+3*panic_threshold; s2ps3=s2+s3;}
    T s1ms2 = s1-s2; if(fabs(s1ms2)<panic_threshold) {s1=s1+1*panic_threshold; s2=s2+2*panic_threshold; s1ms2=s1-s2;}
    
    
    //Now we get a1,a2,a3,a11,a12,a13,a22,a23,a33
    
    
    
    T den = (3*a*a+2*a*(s1+s2+s3)+(s1*s2+s1*s3+s2*s3));if(fabs(den)<panic_threshold) den=den<0?-panic_threshold:panic_threshold;
    
    
       //T s1s1ms2s2 = s1ps2*s1ms2;
    //T/ s1s1ms3s3 = s1ps3*s1ms3;
    //T s2s2ms3s3 = s2ps3*s2ms3;
    
    a1 = -(a+s2)*(a+s3)/den;
    a2 = -(a+s1)*(a+s3)/den;
    a3 = -(a+s1)*(a+s2)/den;
    a11 = -2*a1*(s3+2*a+a1*s3+3*a1*a+s2+a1*s2+a1*s1)/den;
    a22 = -2*a2*(s3+2*a+a2*s3+3*a2*a+s1+a2*s1+a2*s2)/den;
    a33 = -2*a3*(s1+2*a+a3*s1+3*a3*a+s2+a3*s2+a3*s3)/den;
    a12 = -(s3+a+a1*s3+2*a1*a+a1*s1+2*a1*a2*s3+6*a1*a2*a+2*a1*a2*s2+2*a1*a2*s1+a2*s3+2*a2*a+a2*s2)/den;
    a13 = -(s2+a+a1*s2+2*a1*a+a1*s1+2*a1*a3*s2+6*a1*a3*a+2*a1*a3*s3+2*a1*a3*s1+a3*s2+2*a3*a+a3*s3)/den;
    a23 = -(s1+a+a3*s1+2*a3*a+a3*s3+2*a3*a2*s1+6*a3*a2*a+2*a3*a2*s2+2*a3*a2*s3+a2*s1+2*a2*a+a2*s2)/den;
    
    /////////////Hopefully I can get the member function to work at some point
    
    if(triangle==-2){ //This is a signal to print out a and its derivatives
        printf("a  : %10.5e \n",a);
        printf("a1 : %10.5e \n",a1);
        printf("a2 : %10.5e \n",a2);
        printf("a3 : %10.5e \n",a3);
        printf("a11: %10.5e \n",a11);
        printf("a12: %10.5e \n",a12);
        printf("a13: %10.5e \n",a13);
        printf("a22: %10.5e \n",a22);
        printf("a23: %10.5e \n",a23);
        printf("a33: %10.5e \n",a33);
        printf("s1 : %10.5e \n",s1);
        printf("s2 : %10.5e \n",s2);
        printf("s3 : %10.5e \n",s3);    
        printf("mu : %10.5e \n",mu);
        printf("la : %10.5e \n",la);
        printf("c  : %10.5e \n",c);    
        printf("k  : %10.5e \n",k);
    }    
    
    T t362, t363, t364, t365, t366, t367, t368, t369, t370, t371, t372, t373, t374, t375, t376, t377, t378, t379, t380, t381, t382, t383, t384, t385, t386, t387, t388, t389, t390, t391, t392, t393, t394, t395, t396, t397, t398, t399, t400, t401, t402, t403, t404, t405, t406, t407, t408, t409, t410, t411, t412, t413, t414, t415, t416, t417, t418, t419, t420, t421, t422, t423, t424, t425, t426, t427, t428, t429, t430, t431, t432, t433, t434, t435, t436, t437, t438, t439, t440, t441, t442, t443, t444, t445, t446, t447, t448, t449, t450, t451, t452, t453, t454, t455, t456, t457, t458, t459, t460, t461, t462, t463, t464, t465, t466, t467, t468, t469, t470, t471, t472, t473, t474, t475, t476, t477, t478, t479, t480, t481, t482, t483, t484, t485, t486, t487, t488, t489, t490, t491, t492, t493, t494, t495, t496, t497, t498, t499, t500, t501, t502, t503, t504, t505, t506, t507, t508, t509, t510, t511, t512, t513, t514, t515, t516, t517, t518, t519, t520, t521, t522, t523, t524, t525, t526, t527, t528, t529, t530, t531, t532, t533, t534, t535, t536, t537, t538, t539, t540, t541, t542, t543, t544, t545, t546, t547, t548, t549, t550, t551, t552, t553, t554, t555, t556, t557, t558, t559, t560, t561, t562, t563, t564, t565, t566, t567, t568, t569, t570, t571, t572, t573, t574, t575, t576, t577, t578, t579, t580, t581, t582, t583, t584, t585, t586, t587, t588, t589, t590, t591, t592, t593, t594, t595, t596, t597, t598, t599, t600, t601, t602, t603, t604, t605, t606, t607, t608, t609, t610, t611, t612, t613, t614, t615, t616, t617, t618, t619, t620, t621, t622, t623, t624, t625, t626, t627, t628, t629, t630, t631, t632, t633, t634, t635, t636, t637, t638, t639, t640, t641, t642, t643, t644, t645, t646, t647, t648, t649, t650, t651, t652, t653, t654, t655, t656, t657, t658, t659, t660, t661, t662, t663, t664, t665, t666, t667, t668, t669, t670, t671;//, t672, t673, t674;
    
    t362 = a1+1;
    t363 = a+s1;
    t364 = 1/(t363*t363);
    t365 = a+s2;
    t366 = a+s3;
    t367 = t363*t365*t366;
    t368 = log(t367);
    t369 = 1/t363;
    t370 = 1/t365;
    t371 = 1/t366;
    t372 = a*2;
    t373 = a1*a1;
    t374 = 1/(t365*t365);
    t375 = 1/(t365*t365*t365);
    t376 = s1*2;
    t377 = t372+t376;
    t378 = s2*2;
    t379 = t372+t378;
    t380 = s3*2;
    t381 = t372+t380;
    t382 = a1*t363*t365;
    t383 = a1*t363*t366;
    t384 = t362*t365*t366;
    t385 = t382+t383+t384;
    t386 = 1/(t366*t366);
    t387 = 1/(t366*t366*t366);
    t388 = t373*t377;
    t389 = a11*t363*t365;
    t390 = a11*t363*t366;
    t391 = a11*t365*t366;
    t392 = a1*t362*t379;
    t393 = a1*t362*t381;
    t394 = t388+t389+t390+t391+t392+t393;
    t395 = a1*la*t369*t374*t385*t386;
    t396 = a1*mu;
    t397 = t362*t362;
    t398 = 1/(t363*t363*t363);
    t399 = a2+1;
    t400 = mu*t363;
    t401 = la*t368*t369;
    t425 = mu*t369;
    t402 = t400+t401-t425;
    t403 = mu*t365;
    t404 = la*t368*t370;
    t426 = mu*t370;
    t405 = t403+t404-t426;
    t406 = mu*t366;
    t407 = la*t368*t371;
    t427 = mu*t371;
    t408 = t406+t407-t427;
    t409 = a2*a2;
    t410 = a2*t363*t365;
    t411 = a2*t365*t366;
    t412 = t363*t366*t399;
    t413 = t410+t411+t412;
    t414 = t379*t409;
    t415 = a22*t363*t365;
    t416 = a22*t363*t366;
    t417 = a22*t365*t366;
    t418 = a2*t377*t399;
    t419 = a2*t381*t399;
    t420 = t414+t415+t416+t417+t418+t419;
    t421 = a2*la*t364*t370*t386*t413;
    t422 = a2*mu;
    t423 = t399*t399;
    t424 = a3+1;
    t428 = a3*a3;
    t429 = a3*t363*t366;
    t430 = a3*t365*t366;
    t431 = t363*t365*t424;
    t432 = t429+t430+t431;
    t433 = t381*t428;
    t434 = a33*t363*t365;
    t435 = a33*t363*t366;
    t436 = a33*t365*t366;
    t437 = a3*t377*t424;
    t438 = a3*t379*t424;
    t439 = t433+t434+t435+t436+t437+t438;
    t440 = a3*la*t364*t371*t374*t432;
    t441 = a3*mu;
    t442 = t424*t424;
    t443 = mu*t362;
    t444 = mu*t362*t364;
    t445 = la*t364*t370*t371*t385;
    t473 = la*t362*t364*t368;
    t446 = t443+t444+t445-t473;
    t447 = mu*t399;
    t448 = mu*t374*t399;
    t449 = la*t369*t371*t374*t413;
    t498 = la*t368*t374*t399;
    t450 = t447+t448+t449-t498;
    t451 = t362*t366*t399;
    t452 = a12*t363*t365;
    t453 = a12*t363*t366;
    t454 = a12*t365*t366;
    t455 = a1*a2*t363;
    t456 = a1*a2*t365;
    t457 = a1*a2*t366;
    t458 = a1*t363*t399;
    t459 = a2*t362*t365;
    t460 = t451+t452+t453+t454+t455+t456+t457+t458+t459;
    t461 = a1*mu*t374;
    t462 = la*t369*t371*t374*t385;
    t490 = a1*la*t368*t374;
    t463 = t396+t461+t462-t490;
    t464 = a1*mu*t386;
    t465 = la*t369*t370*t385*t386;
    t491 = a1*la*t368*t386;
    t466 = t396+t464+t465-t491;
    t467 = a2*mu*t364;
    t468 = la*t364*t370*t371*t413;
    t514 = a2*la*t364*t368;
    t469 = t422+t467+t468-t514;
    t470 = a2*mu*t386;
    t471 = la*t369*t370*t386*t413;
    t515 = a2*la*t368*t386;
    t472 = t422+t470+t471-t515;
    t474 = mu*t424;
    t475 = mu*t386*t424;
    t476 = la*t369*t370*t386*t432;
    t499 = la*t368*t386*t424;
    t477 = t474+t475+t476-t499;
    t478 = a1*2;
    t479 = t478+2;
    t480 = t362*t365*t424;
    t481 = a13*t363*t365;
    t482 = a13*t363*t366;
    t483 = a13*t365*t366;
    t484 = a1*a3*t363;
    t485 = a1*a3*t365;
    t486 = a1*a3*t366;
    t487 = a1*t363*t424;
    t488 = a3*t362*t366;
    t489 = t480+t481+t482+t483+t484+t485+t486+t487+t488;
    t492 = a3*mu*t364;
    t493 = la*t364*t370*t371*t432;
    t516 = a3*la*t364*t368;
    t494 = t441+t492+t493-t516;
    t495 = a3*mu*t374;
    t496 = la*t369*t371*t374*t432;
    t517 = a3*la*t368*t374;
    t497 = t441+t495+t496-t517;
    t500 = a3*2;
    t501 = t500+2;
    t502 = a2*2;
    t503 = t502+2;
    t504 = t363*t399*t424;
    t505 = a23*t363*t365;
    t506 = a23*t363*t366;
    t507 = a23*t365*t366;
    t508 = a2*a3*t363;
    t509 = a2*a3*t365;
    t510 = a2*a3*t366;
    t511 = a2*t365*t424;
    t512 = a3*t366*t399;
    t513 = t504+t505+t506+t507+t508+t509+t510+t511+t512;
    t518 = s2*s2;
    t519 = s1*s1;
    t520 = t519*t519;
    t521 = a*a;
    t522 = t518*t518;
    t523 = t521*t521;
    t524 = s1*t521;
    t525 = s2*t521;
    t526 = s3*t521;
    t527 = a*s1*s2;
    t528 = a*s1*s3;
    t529 = a*s2*s3;
    t530 = s1*s2*s3;
    t531 = t524+t525+t526+t527+t528+t529+t530+a*t521;
    t532 = log(t531);
    t533 = s3*s3;
    t534 = t520*t521;
    t535 = a*s1*t519*t521*2;
    t536 = t519*t523;
    t537 = t533*t533;
    t538 = mu*t520*t521;
    t539 = a*mu*s1*t519*t521*2;
    t540 = mu*t519*t523;
    t541 = la*t519*t521*t532;
    t542 = a*a1*k*s1*t523*2;
    t543 = a*a1*k*s1*t519*t521*2;
    t544 = a1*k*t519*t523*4;
    t545 = a*la*s1*t521*t532*2;
    t546 = a*a1*la*s1*t519*t532;
    t547 = a*a1*la*s1*t521*t532*2;
    t548 = a1*la*t519*t521*t532*2;
    t549 = t363*t363;
    t550 = a1*s1;
    t583 = la*t532;
    t551 = la+mu-t583;
    t552 = a*s1;
    t553 = a1*t519*2;
    t554 = a*a1*s1*4;
    t555 = a*mu*s2*t521*2;
    t556 = a*mu*s3*t521*2;
    t557 = la*t518*t521;
    t558 = la*t521*t533;
    t559 = mu*t518*t521;
    t560 = mu*t521*t533;
    t561 = a*la*s2*t521*2;
    t562 = a*la*s3*t521*2;
    t563 = a*a2*la*s2*t518;
    t564 = a*a2*la*s2*t521*4;
    t565 = a*a3*la*s3*t533;
    t566 = a*a3*la*s3*t521*4;
    t567 = a*a2*mu*s2*t518;
    t568 = a*a2*mu*s2*t521*2;
    t569 = a*a2*mu*s2*t523*3;
    t570 = a*a3*mu*s3*t533;
    t571 = a*a3*mu*s3*t521*2;
    t572 = a*a3*mu*s3*t523*3;
    t573 = a2*la*t518*t521*4;
    t574 = a3*la*t521*t533*4;
    t575 = a2*mu*t518*t521*2;
    t576 = a*a2*mu*s2*t518*t521*3;
    t577 = a2*mu*t518*t523*6;
    t578 = a3*mu*t521*t533*2;
    t579 = a*a3*mu*s3*t521*t533*3;
    t580 = a3*mu*t523*t533*6;
    t581 = t365*t365;
    t582 = t366*t366;
    t584 = a*mu*s2*t519*t521*2;
    t585 = a*mu*s1*t519*t533*2;
    t586 = mu*s1*s3*t519*t521*4;
    t587 = a*mu*s2*t521*t533*2;
    t588 = a*mu*s1*t523*2;
    t589 = a*la*s1*t519*t521*2;
    t590 = la*t519*t523*5;
    t591 = mu*t519*t523*4;
    t592 = a*mu*s3*t521*t533*2;
    t593 = a*la*s1*t523*3;
    t594 = a*mu*s1*t518*t519*2;
    t595 = mu*s1*s2*t519*t521*4;
    t596 = a*mu*s1*t521*t533*2;
    t597 = a*mu*s3*t519*t521*2;
    t598 = a*mu*s3*t518*t533*2;
    t599 = mu*s2*s3*t521*t533*4;
    t600 = a*la*s1*t368*t518*t521*2;
    t601 = la*s2*s3*t519*t521*4;
    t602 = mu*s2*s3*t519*t521*4;
    t603 = a1*k*t518*t523*t533*4;
    t604 = la*t368*t518*t521*t533*4;
    t605 = a1*k*s2*s3*t521*t523*4;
    t606 = la*s2*s3*t368*t523*4;
    t607 = a*la*s1*s2*s3*t519;
    t608 = a*mu*s1*s2*s3*t519*2;
    t609 = a*a1*k*s1*t518*t521*t533*8;
    t610 = a1*la*s2*s3*t368*t523*4;
    t611 = a*a1*k*s1*s2*s3*t523*8;
    t612 = a1*k*t518*t519*t521*t533*4;
    t613 = a1*k*s2*s3*t519*t523*4;
    t614 = a*la*s1*t368*t518*t533*2;
    t615 = a1*la*t368*t518*t521*t533*2;
    t616 = a1*la*s2*s3*t368*t519*t521*2;
    t617 = a*a1*la*s1*s2*s3*t368*t521*4;
    t618 = la*t518*t519*t521*3;
    t619 = la*t519*t521*t533*3;
    t620 = mu*t518*t519*t521*4;
    t621 = mu*t519*t521*t533*4;
    t622 = a*la*s2*t368*t523*2;
    t623 = a*la*s3*t368*t523*2;
    t624 = la*s1*s2*t523*5;
    t625 = la*s1*s3*t523*5;
    t626 = mu*s1*s2*t523*4;
    t627 = mu*s1*s3*t523*4;
    t628 = a*la*s2*t368*t518*t521*2;
    t629 = la*t368*t518*t523*4;
    t630 = a*la*s3*t368*t521*t533*2;
    t631 = la*t368*t523*t533*4;
    t632 = a*la*s2*t368*t521*t533*2;
    t633 = a*la*s3*t368*t518*t521*2;
    t634 = a*la*s2*t519*t533;
    t635 = a*la*s3*t518*t519;
    t636 = a*mu*s2*t519*t533*2;
    t637 = a*mu*s3*t518*t519*2;
    t638 = a2*la*t519*t521*t533*4;
    t639 = a3*la*t518*t519*t521*4;
    t640 = a2*mu*t519*t521*t533*2;
    t641 = a3*mu*t518*t519*t521*2;
    t642 = a2*mu*t519*t523*t533*6;
    t643 = a3*mu*t518*t519*t523*6;
    t644 = a2*la*s1*s3*t523*12;
    t645 = a3*la*s1*s2*t523*12;
    t646 = a2*mu*s1*s3*t523*4;
    t647 = a3*mu*s1*s2*t523*4;
    t648 = a2*mu*s1*s3*t521*t523*6;
    t649 = a3*mu*s1*s2*t521*t523*6;
    t650 = a*a2*mu*s2*t519*t521*t533*12;
    t651 = a*a3*mu*s3*t518*t519*t521*12;
    t652 = a*a2*la*s1*s2*s3*t521*14;
    t653 = a*a3*la*s1*s2*s3*t521*14;
    t654 = a*a2*mu*s1*s2*s3*t521*4;
    t655 = a*a2*mu*s1*s2*s3*t523*12;
    t656 = a*a3*mu*s1*s2*s3*t521*4;
    t657 = a*a3*mu*s1*s2*s3*t523*12;
    t658 = a*la*s1*s2*s3*t368*t533*2;
    t659 = a*la*s1*s2*s3*t368*t518*2;
    t660 = a2*mu*t518*t519*t521*t533*6;
    t661 = a3*mu*t518*t519*t521*t533*6;
    t662 = a*a2*la*s2*t519*t533*2;
    t663 = a2*la*s1*s3*t518*t521*4;
    t664 = a*a3*la*s3*t518*t519*2;
    t665 = a3*la*s1*s2*t521*t533*4;
    t666 = a2*mu*s1*s3*t518*t521*2;
    t667 = a2*mu*s1*s3*t518*t523*6;
    t668 = a3*mu*s1*s2*t521*t533*2;
    t669 = a3*mu*s1*s2*t523*t533*6;
    t670 = la*s1*s2*t368*t521*t533*4;
    t671 = la*s1*s3*t368*t518*t521*4;
    dP_dF.x1111 = a*(t395-a11*mu-a11*mu*t374+mu*t373*t375*2+a11*la*t368*t374-la*t368*t373*t375*2-la*t369*t371*t374*t394+la*t362*t364*t371*t374*t385+a1*la*t369*t371*t375*t385*3)+a*(t395-a11*mu-a11*mu*t386+mu*t373*t387*2+a11*la*t368*t386-la*t368*t373*t387*2-la*t369*t370*t386*t394+la*t362*t364*t370*t385*t386+a1*la*t369*t370*t385*t387*3)-a1*t446*2-a1*t463*2-a1*t466*2-a11*t402-a11*t405-a11*t408+k*t373*2+mu*(t373*4+t397*2+a11*t377+a11*t379+a11*t381)/2+a*a11*k*2-a*a11*mu-a*a11*mu*t364+a*mu*t397*t398*2+a*a11*la*t364*t368-a*la*t368*t397*t398*2-mu*t369*t370*t371*t394+la*t364*t374*(t385*t385)*t386+la*t368*t369*t370*t371*t394+mu*t362*t364*t370*t371*t385-a*la*t364*t370*t371*t394+a1*mu*t369*t371*t374*t385+a1*mu*t369*t370*t385*t386+a*a1*la*t364*t371*t374*t385+a*a1*la*t364*t370*t385*t386+a*la*t362*t370*t371*t385*t398*3-a1*la*t368*t369*t371*t374*t385-a1*la*t368*t369*t370*t385*t386-la*t362*t364*t368*t370*t371*t385;
    dP_dF.x2222 = a*(t421-a22*mu-a22*mu*t364+mu*t398*t409*2+a22*la*t364*t368-la*t368*t398*t409*2-la*t364*t370*t371*t420+la*t364*t371*t374*t399*t413+a2*la*t370*t371*t398*t413*3)+a*(t421-a22*mu-a22*mu*t386+mu*t387*t409*2+a22*la*t368*t386-la*t368*t387*t409*2-la*t369*t370*t386*t420+la*t369*t374*t386*t399*t413+a2*la*t369*t370*t387*t413*3)-a2*t450*2-a2*t469*2-a2*t472*2-a22*t402-a22*t405-a22*t408+k*t409*2+mu*(t409*4+t423*2+a22*t377+a22*t379+a22*t381)/2+a*a22*k*2-a*a22*mu-a*a22*mu*t374+a*mu*t375*t423*2+a*a22*la*t368*t374-a*la*t368*t375*t423*2-mu*t369*t370*t371*t420+la*t364*t374*t386*(t413*t413)+la*t368*t369*t370*t371*t420+mu*t369*t371*t374*t399*t413-a*la*t369*t371*t374*t420+a2*mu*t364*t370*t371*t413+a2*mu*t369*t370*t386*t413+a*a2*la*t364*t371*t374*t413+a*a2*la*t369*t374*t386*t413+a*la*t369*t371*t375*t399*t413*3-a2*la*t364*t368*t370*t371*t413-a2*la*t368*t369*t370*t386*t413-la*t368*t369*t371*t374*t399*t413;
    dP_dF.x3333 = a*(t440-a33*mu-a33*mu*t374+mu*t375*t428*2+a33*la*t368*t374-la*t368*t375*t428*2-la*t369*t371*t374*t439+la*t369*t374*t386*t424*t432+a3*la*t369*t371*t375*t432*3)+a*(t440-a33*mu-a33*mu*t364+mu*t398*t428*2+a33*la*t364*t368-la*t368*t398*t428*2-la*t364*t370*t371*t439+la*t364*t370*t386*t424*t432+a3*la*t370*t371*t398*t432*3)-a3*t477*2-a3*t494*2-a3*t497*2-a33*t402-a33*t405-a33*t408+k*t428*2+mu*(t428*4+t442*2+a33*t377+a33*t379+a33*t381)/2+a*a33*k*2-a*a33*mu-a*a33*mu*t386+a*mu*t387*t442*2+a*a33*la*t368*t386-a*la*t368*t387*t442*2-mu*t369*t370*t371*t439+la*t364*t374*t386*(t432*t432)+la*t368*t369*t370*t371*t439+mu*t369*t370*t386*t424*t432-a*la*t369*t370*t386*t439+a3*mu*t364*t370*t371*t432+a3*mu*t369*t371*t374*t432+a*a3*la*t364*t370*t386*t432+a*a3*la*t369*t374*t386*t432+a*la*t369*t370*t387*t424*t432*3-a3*la*t364*t368*t370*t371*t432-a3*la*t368*t369*t371*t374*t432-la*t368*t369*t370*t386*t424*t432;
    dP_dF.x2211 = a*(-a12*mu-a12*mu*t386+a1*a2*mu*t387*2+a12*la*t368*t386-a1*a2*la*t368*t387*2-la*t369*t370*t386*t460+la*t369*t374*t385*t386*t399+a1*la*t369*t370*t387*t413+a2*la*t364*t370*t385*t386+a2*la*t369*t370*t385*t387*2)-a1*t450-a1*t469-a1*t472-a12*t402-a12*t405-a12*t408-a2*t446-a2*t463-a2*t466+a*(-a12*mu-a12*mu*t374+a12*la*t368*t374+a1*mu*t375*t399*2-a1*la*t368*t375*t399*2-la*t369*t371*t374*t460+la*t369*t371*t375*t385*t399*2+a1*la*t369*t371*t375*t413+a2*la*t364*t371*t374*t385+a2*la*t369*t374*t385*t386)+mu*(a1*a2*2+a12*t377+a12*t379+a12*t381+a1*t503+a2*t479)/2+a*a12*k*2+a1*a2*k*2-a*a12*mu-a*a12*mu*t364+a*a12*la*t364*t368+a*a2*mu*t362*t398*2-mu*t369*t370*t371*t460+la*t364*t374*t385*t386*t413+la*t368*t369*t370*t371*t460+mu*t369*t371*t374*t385*t399-a*a2*la*t362*t368*t398*2-a*la*t364*t370*t371*t460+a2*mu*t364*t370*t371*t385+a2*mu*t369*t370*t385*t386+a*a2*la*t364*t370*t385*t386+a*a2*la*t370*t371*t385*t398*2+a*la*t364*t371*t374*t385*t399+a*la*t362*t370*t371*t398*t413-a2*la*t364*t368*t370*t371*t385-a2*la*t368*t369*t370*t385*t386-la*t368*t369*t371*t374*t385*t399;
    dP_dF.x3311 = a*(-a13*mu-a13*mu*t374+a1*a3*mu*t375*2+a13*la*t368*t374-a1*a3*la*t368*t375*2-la*t369*t371*t374*t489+la*t369*t374*t385*t386*t424+a1*la*t369*t371*t375*t432+a3*la*t364*t371*t374*t385+a3*la*t369*t371*t375*t385*2)-a1*t477-a1*t494-a1*t497-a13*t402-a13*t405-a13*t408-a3*t446-a3*t463-a3*t466+a*(-a13*mu-a13*mu*t386+a13*la*t368*t386+a1*mu*t387*t424*2-a1*la*t368*t387*t424*2-la*t369*t370*t386*t489+la*t369*t370*t385*t387*t424*2+a1*la*t369*t370*t387*t432+a3*la*t364*t370*t385*t386+a3*la*t369*t374*t385*t386)+mu*(a1*a3*2+a13*t377+a13*t379+a13*t381+a1*t501+a3*t479)/2+a*a13*k*2+a1*a3*k*2-a*a13*mu-a*a13*mu*t364+a*a13*la*t364*t368+a*a3*mu*t362*t398*2-mu*t369*t370*t371*t489+la*t364*t374*t385*t386*t432+la*t368*t369*t370*t371*t489+mu*t369*t370*t385*t386*t424-a*a3*la*t362*t368*t398*2-a*la*t364*t370*t371*t489+a3*mu*t364*t370*t371*t385+a3*mu*t369*t371*t374*t385+a*a3*la*t364*t371*t374*t385+a*a3*la*t370*t371*t385*t398*2+a*la*t364*t370*t385*t386*t424+a*la*t362*t370*t371*t398*t432-a3*la*t364*t368*t370*t371*t385-a3*la*t368*t369*t371*t374*t385-la*t368*t369*t370*t385*t386*t424;
    dP_dF.x3322 = a*(-a23*mu-a23*mu*t364+a2*a3*mu*t398*2+a23*la*t364*t368-a2*a3*la*t368*t398*2-la*t364*t370*t371*t513+la*t364*t370*t386*t413*t424+a2*la*t370*t371*t398*t432+a3*la*t364*t371*t374*t413+a3*la*t370*t371*t398*t413*2)-a2*t477-a2*t494-a2*t497-a23*t402-a23*t405-a23*t408-a3*t450-a3*t469-a3*t472+a*(-a23*mu-a23*mu*t386+a23*la*t368*t386+a2*mu*t387*t424*2-a2*la*t368*t387*t424*2-la*t369*t370*t386*t513+la*t369*t370*t387*t413*t424*2+a2*la*t369*t370*t387*t432+a3*la*t364*t370*t386*t413+a3*la*t369*t374*t386*t413)+mu*(a2*a3*2+a23*t377+a23*t379+a23*t381+a2*t501+a3*t503)/2+a*a23*k*2+a2*a3*k*2-a*a23*mu-a*a23*mu*t374+a*a23*la*t368*t374+a*a3*mu*t375*t399*2-mu*t369*t370*t371*t513+la*t364*t374*t386*t413*t432+la*t368*t369*t370*t371*t513+mu*t369*t370*t386*t413*t424-a*a3*la*t368*t375*t399*2-a*la*t369*t371*t374*t513+a3*mu*t364*t370*t371*t413+a3*mu*t369*t371*t374*t413+a*a3*la*t364*t371*t374*t413+a*a3*la*t369*t371*t375*t413*2+a*la*t369*t371*t375*t399*t432+a*la*t369*t374*t386*t413*t424-a3*la*t364*t368*t370*t371*t413-a3*la*t368*t369*t371*t374*t413-la*t368*t369*t370*t386*t413*t424;
    dP_dF.x2121 = (t538+t539+t540+t541+t542+t543+t544+t545+t546+t547+t548+t555+t557+t559+t561+t563+t564+t567+t568+t569+t573+t575+t576+t577+t584+t594+t595-la*t519*t521+mu*t518*t520-mu*t519*t521-mu*t518*t523-mu*t519*t522-mu*t521*t522-a2*k*t518*t523*4-a*la*s1*t521*2-a*mu*s1*t521*2-a*mu*s1*t522*2+a*mu*s2*t520*2-a1*la*t519*t521*4-a1*mu*t519*t521*2-a1*mu*t519*t523*6-la*t518*t521*t532-a*a2*k*s2*t523*2-a*a1*la*s1*t518-a*a1*la*s1*t519-a*a1*la*s1*t521*4+a*a2*la*s1*t518*2-a*a1*la*s2*t519*2+a*a2*la*s2*t519-a*a1*mu*s1*t518-a*a1*mu*s1*t519-a*a1*mu*s1*t521*2-a*a1*mu*s1*t523*3+a*a2*mu*s2*t519+a1*k*s1*s2*t523*4-a2*k*s1*s2*t523*4-a1*la*s1*s2*t521*4+a2*la*s1*s2*t521*4+a1*k*t518*t519*t521*4-a2*k*t518*t519*t521*4-a*la*s2*t521*t532*2-a1*mu*s1*s2*t521*2-a1*mu*s1*s2*t523*6+a2*mu*s1*s2*t521*2+a2*mu*s1*s2*t523*6-a*mu*s1*t518*t521*2-a*mu*s2*t518*t519*2-a*mu*s2*t518*t521*2-a2*la*t518*t521*t532*2-a1*mu*t518*t519*t521*6+a2*mu*t518*t519*t521*6-mu*s1*s2*t518*t521*4+a*a1*k*s1*t518*t519*2+a*a1*k*s1*t518*t521*2-a*a2*k*s1*t518*t521*8+a*a1*k*s2*t519*t521*8-a*a2*k*s2*t518*t519*2-a*a2*k*s2*t518*t521*2-a*a2*k*s2*t519*t521*2+a*a1*la*s1*t518*t532-a*a2*la*s2*t518*t532-a*a2*la*s2*t519*t532-a*a2*la*s2*t521*t532*2-a*a1*mu*s1*t518*t519*3-a*a1*mu*s1*t518*t521*3-a*a1*mu*s1*t519*t521*3+a*a2*mu*s1*t518*t521*12-a*a1*mu*s2*t519*t521*12+a*a2*mu*s2*t518*t519*3+a*a2*mu*s2*t519*t521*3+a1*k*s1*s2*t519*t521*4-a2*k*s1*s2*t518*t521*4+a1*la*s1*s2*t521*t532*2-a2*la*s1*s2*t521*t532*2-a1*mu*s1*s2*t519*t521*6+a2*mu*s1*s2*t518*t521*6)/(t534+t535+t536+t518*t520-t518*t523-t519*t522-t521*t522-a*s1*t522*2+a*s2*t520*2+a*s1*t518*t519*2-a*s1*t518*t521*2-a*s2*t518*t519*2-a*s2*t518*t521*2+a*s2*t519*t521*2-s1*s2*t518*t521*4+s1*s2*t519*t521*4)-(t364*t374*(a*la*t363*t365*t371*(t552+t553+t554-a*s2-a2*t518*2-a*a2*s2*4+a1*s1*s2*2-a2*s1*s2*2)+a*t386*t549*t551*t581*(t550-a2*s2)))/((s1+s2)*(s1-s2));
    dP_dF.x3131 = (t538+t539+t540+t541+t542+t543+t544+t545+t546+t547+t548+t556+t558+t560+t562+t565+t566+t570+t571+t572+t574+t578+t579+t580+t585+t586+t597-la*t519*t521-mu*t519*t521-mu*t519*t537+mu*t520*t533-mu*t523*t533-mu*t521*t537-a3*k*t523*t533*4-a*la*s1*t521*2-a*mu*s1*t521*2-a*mu*s1*t537*2+a*mu*s3*t520*2-a1*la*t519*t521*4-a1*mu*t519*t521*2-a1*mu*t519*t523*6-la*t521*t532*t533-a*a3*k*s3*t523*2-a*a1*la*s1*t519-a*a1*la*s1*t521*4-a*a1*la*s1*t533-a*a1*la*s3*t519*2+a*a3*la*s1*t533*2+a*a3*la*s3*t519-a*a1*mu*s1*t519-a*a1*mu*s1*t521*2-a*a1*mu*s1*t523*3-a*a1*mu*s1*t533+a*a3*mu*s3*t519+a1*k*s1*s3*t523*4-a3*k*s1*s3*t523*4-a1*la*s1*s3*t521*4+a3*la*s1*s3*t521*4+a1*k*t519*t521*t533*4-a3*k*t519*t521*t533*4-a*la*s3*t521*t532*2-a1*mu*s1*s3*t521*2-a1*mu*s1*s3*t523*6+a3*mu*s1*s3*t521*2+a3*mu*s1*s3*t523*6-a*mu*s1*t521*t533*2-a*mu*s3*t519*t533*2-a*mu*s3*t521*t533*2-a3*la*t521*t532*t533*2-a1*mu*t519*t521*t533*6+a3*mu*t519*t521*t533*6-mu*s1*s3*t521*t533*4+a*a1*k*s1*t519*t533*2+a*a1*k*s1*t521*t533*2+a*a1*k*s3*t519*t521*8-a*a3*k*s1*t521*t533*8-a*a3*k*s3*t519*t521*2-a*a3*k*s3*t519*t533*2-a*a3*k*s3*t521*t533*2+a*a1*la*s1*t532*t533-a*a3*la*s3*t519*t532-a*a3*la*s3*t521*t532*2-a*a3*la*s3*t532*t533-a*a1*mu*s1*t519*t521*3-a*a1*mu*s1*t519*t533*3-a*a1*mu*s1*t521*t533*3-a*a1*mu*s3*t519*t521*12+a*a3*mu*s1*t521*t533*12+a*a3*mu*s3*t519*t521*3+a*a3*mu*s3*t519*t533*3+a1*k*s1*s3*t519*t521*4-a3*k*s1*s3*t521*t533*4+a1*la*s1*s3*t521*t532*2-a3*la*s1*s3*t521*t532*2-a1*mu*s1*s3*t519*t521*6+a3*mu*s1*s3*t521*t533*6)/(t534+t535+t536-t519*t537+t520*t533-t523*t533-t521*t537-a*s1*t537*2+a*s3*t520*2+a*s1*t519*t533*2-a*s1*t521*t533*2+a*s3*t519*t521*2-a*s3*t519*t533*2-a*s3*t521*t533*2+s1*s3*t519*t521*4-s1*s3*t521*t533*4)-(t364*t386*(a*la*t363*t366*t370*(t552+t553+t554-a*s3-a3*t533*2-a*a3*s3*4+a1*s1*s3*2-a3*s1*s3*2)+a*t374*t549*t551*t582*(t550-a3*s3)))/((s1+s3)*(s1-s3));
    dP_dF.x3232 = -(t555-t556+t557-t558+t559-t560+t561-t562+t563+t564-t565-t566+t567+t568+t569-t570-t571-t572+t573-t574+t575+t576+t577-t578-t579-t580+t587+t592+t598+t599-mu*t518*t523-mu*t521*t522+mu*t518*t537-mu*t522*t533+mu*t523*t533+mu*t521*t537-a2*k*t518*t523*4+a3*k*t523*t533*4+a*mu*s2*t537*2-a*mu*s3*t522*2-la*t518*t521*t532+la*t521*t532*t533-a*a2*k*s2*t523*2+a*a3*k*s3*t523*2+a*a2*la*s2*t533+a*a2*la*s3*t518*2-a*a3*la*s2*t533*2-a*a3*la*s3*t518+a*a2*mu*s2*t533-a*a3*mu*s3*t518-a2*k*s2*s3*t523*4+a3*k*s2*s3*t523*4+a2*la*s2*s3*t521*4-a3*la*s2*s3*t521*4-a2*k*t518*t521*t533*4+a3*k*t518*t521*t533*4-a*la*s2*t521*t532*2+a*la*s3*t521*t532*2+a2*mu*s2*s3*t521*2+a2*mu*s2*s3*t523*6-a3*mu*s2*s3*t521*2-a3*mu*s2*s3*t523*6-a*mu*s2*t518*t521*2-a*mu*s2*t518*t533*2-a*mu*s3*t518*t521*2-a2*la*t518*t521*t532*2+a3*la*t521*t532*t533*2+a2*mu*t518*t521*t533*6-a3*mu*t518*t521*t533*6-mu*s2*s3*t518*t521*4-a*a2*k*s2*t518*t521*2-a*a2*k*s2*t518*t533*2-a*a2*k*s2*t521*t533*2-a*a2*k*s3*t518*t521*8+a*a3*k*s2*t521*t533*8+a*a3*k*s3*t518*t521*2+a*a3*k*s3*t518*t533*2+a*a3*k*s3*t521*t533*2-a*a2*la*s2*t518*t532-a*a2*la*s2*t521*t532*2-a*a2*la*s2*t532*t533+a*a3*la*s3*t518*t532+a*a3*la*s3*t521*t532*2+a*a3*la*s3*t532*t533+a*a2*mu*s2*t518*t533*3+a*a2*mu*s2*t521*t533*3+a*a2*mu*s3*t518*t521*12-a*a3*mu*s2*t521*t533*12-a*a3*mu*s3*t518*t521*3-a*a3*mu*s3*t518*t533*3-a2*k*s2*s3*t518*t521*4+a3*k*s2*s3*t521*t533*4-a2*la*s2*s3*t521*t532*2+a3*la*s2*s3*t521*t532*2+a2*mu*s2*s3*t518*t521*6-a3*mu*s2*s3*t521*t533*6)/(t518*t523+t521*t522-t518*t537+t522*t533-t523*t533-t521*t537-a*s2*t537*2+a*s3*t522*2+a*s2*t518*t521*2+a*s2*t518*t533*2-a*s2*t521*t533*2+a*s3*t518*t521*2-a*s3*t518*t533*2-a*s3*t521*t533*2+s2*s3*t518*t521*4-s2*s3*t521*t533*4)-(t374*t386*(a*la*t365*t366*t369*(a*s2-a*s3+a2*t518*2-a3*t533*2+a*a2*s2*4-a*a3*s3*4+a2*s2*s3*2-a3*s2*s3*2)+a*t364*t551*t581*t582*(a2*s2-a3*s3)))/((s2+s3)*(s2-s3));
    dP_dF.x2112 = -(t364*t374*t386*(t539+t584+t585+t586-t587+t588+t589+t590+t591+t593+t596+t600+t601+t602+t603+t604+t605+t606+t607+t608+t609+t610+t611+t612+t613+t614+t615+t616+t617+t619+t621+t622+t625+t627+t628+t629+t632+t634+t636+t638+t640+t642+t644+t646+t648+t650+t652+t654+t655+t659+t660+t662+t663+t666+t667+t671-la*t518*t523*5-mu*t518*t523*4-a*la*s2*t523*3-a*mu*s2*t523*2-a1*la*t518*t523*12+a2*la*t519*t523*12-a1*mu*t518*t523*4+a2*mu*t519*t523*4-la*s2*s3*t523*5-mu*s2*s3*t523*4-la*t368*t519*t523*4-la*t518*t521*t533*3-mu*t518*t521*t533*4-a*a1*la*s2*t523*9+a*a2*la*s1*t523*9-a*a1*mu*s2*t523*3+a*a2*mu*s1*t523*3-a1*la*s1*s2*t523*12+a2*la*s1*s2*t523*12-a1*la*s2*s3*t523*12+a1*k*t518*t519*t523*4+a1*k*t518*t521*t523*4-a2*k*t518*t519*t523*4-a2*k*t519*t521*t523*4-a2*k*t519*t523*t533*4-a*la*s1*t368*t523*2-a*la*s1*t518*t521*3-a*la*s1*t518*t533+a*la*s1*t519*t533+a*la*s1*t521*t533*2-a*la*s2*t518*t521*2+a*la*s2*t519*t521*3-a*la*s2*t518*t533-a*la*s2*t521*t533*2-a*la*s3*t518*t521*8+a*la*s3*t519*t521*8-a1*mu*s1*s2*t523*4+a2*mu*s1*s2*t523*4-a1*mu*s2*s3*t523*4-a*mu*s1*t518*t521*2-a*mu*s1*t518*t533*2-a*mu*s2*t518*t521*2-a*mu*s2*t518*t533*2-a*mu*s3*t518*t521*8+a*mu*s3*t519*t521*8+a1*la*t368*t518*t523*4-a2*la*t368*t519*t523*4-a1*la*t518*t519*t521*4-a1*la*t518*t521*t533*4+a2*la*t518*t519*t521*4-a1*mu*t518*t519*t521*2-a1*mu*t518*t519*t523*6-a1*mu*t518*t521*t523*6-a1*mu*t518*t521*t533*2-a1*mu*t518*t523*t533*6+a2*mu*t518*t519*t521*2+a2*mu*t518*t519*t523*6+a2*mu*t519*t521*t523*6-la*s1*s3*t368*t523*4-la*s1*s2*t518*t521+la*s1*s2*t519*t521-la*s1*s3*t518*t521*4+la*s1*s3*t519*t521*3-la*s2*s3*t518*t521*3-mu*s1*s2*t518*t521+mu*s1*s2*t519*t521-mu*s1*s2*t518*t533+mu*s1*s2*t519*t533-mu*s1*s3*t518*t521*4-mu*s2*s3*t518*t521*4-la*t368*t519*t521*t533*4+la*s1*s2*t368*t518*t521-la*s1*s2*t368*t519*t521+la*s1*s2*t368*t518*t533-la*s1*s2*t368*t519*t533-la*s1*s3*t368*t519*t521*4+la*s2*s3*t368*t518*t521*4-la*s2*s3*t368*t519*t521*4+a*a1*k*s1*t518*t523*8+a*a1*k*s2*t518*t523*2-a*a2*k*s1*t518*t523*2+a*a1*k*s2*t519*t523*2-a*a2*k*s1*t519*t523*2+a*a1*k*s2*t521*t523*2-a*a2*k*s1*t521*t523*2+a*a1*k*s2*t523*t533*2-a*a2*k*s1*t523*t533*2+a*a1*k*s3*t518*t523*8-a*a2*k*s2*t519*t523*8-a*a2*k*s3*t519*t523*8+a*a1*la*s2*t368*t523*3-a*a2*la*s1*t368*t523*3-a*a1*la*s1*t518*t521*14-a*a1*la*s1*t518*t533*2-a*a1*la*s2*t518*t519+a*a2*la*s1*t518*t519-a*a1*la*s2*t518*t521*4+a*a2*la*s1*t518*t521*4-a*a1*la*s2*t519*t521*4+a*a2*la*s1*t519*t521*4-a*a1*la*s2*t518*t533+a*a2*la*s1*t518*t533-a*a1*la*s2*t519*t533+a*a2*la*s1*t519*t533-a*a1*la*s2*t521*t533*4+a*a2*la*s1*t521*t533*4-a*a1*la*s3*t518*t519*2-a*a1*la*s3*t518*t521*14+a*a2*la*s2*t519*t521*14+a*a2*la*s3*t518*t519*2+a*a2*la*s3*t519*t521*14-a*a1*mu*s1*t518*t521*4-a*a1*mu*s1*t518*t523*12-a*a1*mu*s2*t518*t519+a*a2*mu*s1*t518*t519-a*a1*mu*s2*t518*t521*2+a*a2*mu*s1*t518*t521*2-a*a1*mu*s2*t519*t521*2+a*a2*mu*s1*t519*t521*2-a*a1*mu*s2*t518*t523*3+a*a2*mu*s1*t518*t523*3-a*a1*mu*s2*t519*t523*3+a*a2*mu*s1*t519*t523*3-a*a1*mu*s2*t521*t523*3+a*a2*mu*s1*t521*t523*3-a*a1*mu*s2*t518*t533+a*a2*mu*s1*t518*t533-a*a1*mu*s2*t519*t533+a*a2*mu*s1*t519*t533-a*a1*mu*s2*t521*t533*2+a*a2*mu*s1*t521*t533*2-a*a1*mu*s2*t523*t533*3+a*a2*mu*s1*t523*t533*3-a*a1*mu*s3*t518*t521*4+a*a2*mu*s2*t519*t521*4-a*a1*mu*s3*t518*t523*12+a*a2*mu*s2*t519*t523*12+a*a2*mu*s3*t519*t521*4+a*a2*mu*s3*t519*t523*12+a1*k*s1*s2*t518*t523*4+a1*k*s1*s2*t521*t523*4+a1*k*s1*s2*t523*t533*4+a1*k*s1*s3*t518*t523*16-a2*k*s1*s2*t519*t523*4-a2*k*s1*s2*t521*t523*4-a2*k*s1*s2*t523*t533*4+a1*k*s2*s3*t518*t523*4-a2*k*s1*s3*t518*t523*4-a2*k*s1*s3*t519*t523*4-a2*k*s1*s3*t521*t523*4-a2*k*s2*s3*t519*t523*16-a*la*s1*s2*s3*t518-a*mu*s1*s2*s3*t518*2+a1*la*s1*s2*t368*t523*4-a2*la*s1*s2*t368*t523*4-a1*la*s1*s2*t518*t521*4-a1*la*s1*s2*t521*t533*4-a2*la*s1*s3*t368*t523*4-a1*la*s1*s3*t518*t521*12+a2*la*s1*s2*t519*t521*4+a2*la*s1*s2*t521*t533*4-a1*la*s2*s3*t518*t521*4-a1*la*s2*s3*t519*t521*4+a2*la*s1*s3*t519*t521*4+a2*la*s2*s3*t519*t521*12-a2*k*t518*t519*t521*t533*4-a*la*s1*t368*t519*t521*2-a*la*s1*t368*t519*t533*2-a*la*s1*t368*t521*t533*2-a*la*s2*t368*t519*t521*2+a*la*s2*t368*t518*t533*2-a*la*s2*t368*t519*t533*2+a*la*s3*t368*t518*t521*8-a*la*s3*t368*t519*t521*8-a1*mu*s1*s2*t518*t521*2-a1*mu*s1*s2*t518*t523*6-a1*mu*s1*s2*t521*t523*6-a1*mu*s1*s2*t521*t533*2-a1*mu*s1*s2*t523*t533*6+a2*mu*s1*s2*t519*t521*2-a1*mu*s1*s3*t518*t523*24+a2*mu*s1*s2*t519*t523*6+a2*mu*s1*s2*t521*t523*6+a2*mu*s1*s2*t521*t533*2+a2*mu*s1*s2*t523*t533*6-a1*mu*s2*s3*t518*t521*2-a1*mu*s2*s3*t519*t521*2+a2*mu*s1*s3*t519*t521*2-a1*mu*s2*s3*t518*t523*6-a1*mu*s2*s3*t519*t523*6+a2*mu*s1*s3*t519*t523*6-a1*mu*s2*s3*t521*t523*6+a2*mu*s2*s3*t519*t523*24+a1*la*t368*t518*t519*t521*2-a2*la*t368*t518*t519*t521*2-a2*la*t368*t519*t521*t533*2-a1*mu*t518*t519*t521*t533*6-a*a2*k*s1*s2*s3*t523*8-a*a1*la*s1*s2*s3*t518*2-a*a1*la*s1*s2*s3*t521*14+a*a2*la*s1*s2*s3*t519*2+a*a1*k*s2*t518*t519*t521*2-a*a2*k*s1*t518*t519*t521*2+a*a1*k*s2*t518*t519*t533*2-a*a2*k*s1*t518*t519*t533*2+a*a1*k*s2*t518*t521*t533*2-a*a2*k*s1*t518*t521*t533*2+a*a1*k*s2*t519*t521*t533*2-a*a2*k*s1*t519*t521*t533*2+a*a1*k*s3*t518*t519*t521*8-a*a2*k*s2*t519*t521*t533*8-a*a2*k*s3*t518*t519*t521*8-a*a1*mu*s1*s2*s3*t521*4-a*a1*mu*s1*s2*s3*t523*12+a*a1*la*s1*t368*t518*t521*4+a*a1*la*s2*t368*t518*t519-a*a2*la*s1*t368*t518*t519+a*a1*la*s2*t368*t518*t521*2-a*a2*la*s1*t368*t518*t521*2+a*a1*la*s2*t368*t519*t521*2-a*a2*la*s1*t368*t519*t521*2+a*a1*la*s2*t368*t518*t533-a*a2*la*s1*t368*t518*t533+a*a1*la*s2*t368*t519*t533-a*a2*la*s1*t368*t519*t533+a*a1*la*s2*t368*t521*t533*2-a*a2*la*s1*t368*t521*t533*2+a*a1*la*s3*t368*t518*t521*4-a*a2*la*s2*t368*t519*t521*4-a*a2*la*s3*t368*t519*t521*4-a*a1*mu*s1*t518*t521*t533*12-a*a1*mu*s2*t518*t519*t521*3+a*a2*mu*s1*t518*t519*t521*3-a*a1*mu*s2*t518*t519*t533*3+a*a2*mu*s1*t518*t519*t533*3-a*a1*mu*s2*t518*t521*t533*3+a*a2*mu*s1*t518*t521*t533*3-a*a1*mu*s2*t519*t521*t533*3+a*a2*mu*s1*t519*t521*t533*3-a*a1*mu*s3*t518*t519*t521*12+a*a2*mu*s3*t518*t519*t521*12+a1*k*s1*s2*t518*t521*t533*4-a2*k*s1*s2*t519*t521*t533*4+a1*k*s2*s3*t518*t519*t521*4-a2*k*s1*s3*t518*t519*t521*4-a*la*s1*s2*s3*t368*t519*2+a1*la*s1*s2*t368*t518*t521*2+a1*la*s1*s2*t368*t521*t533*2-a2*la*s1*s2*t368*t519*t521*2-a2*la*s1*s2*t368*t521*t533*2+a1*la*s2*s3*t368*t518*t521*2-a2*la*s1*s3*t368*t518*t521*2-a2*la*s1*s3*t368*t519*t521*2-a1*mu*s1*s2*t518*t521*t533*6+a2*mu*s1*s2*t519*t521*t533*6-a1*mu*s2*s3*t518*t519*t521*6+a2*mu*s1*s3*t518*t519*t521*6+a*a1*k*s1*s2*s3*t518*t521*8-a*a2*k*s1*s2*s3*t519*t521*8-a*a2*la*s1*s2*s3*t368*t521*4-a*a1*mu*s1*s2*s3*t518*t521*12+a*a2*mu*s1*s2*s3*t519*t521*12))/(t518-t519);
    

    dP_dF.x3113 = (t364*t374*t386*(t539+t588+t589+t590+t591-t592+t593+t594+t595-t596+t597-t598-t599-t600+t601+t602+t603+t604+t605+t606+t607+t608+t609+t610+t611+t612+t613+t614+t615+t616+t617+t618+t620+t623+t624+t626+t630+t631+t633+t635+t637+t639+t641+t643+t645+t647+t649+t651+t653+t656+t657+t658+t661+t664+t665+t668+t669+t670-la*t523*t533*5-mu*t523*t533*4-a*la*s3*t523*3-a*mu*s3*t523*2-a1*la*t523*t533*12+a3*la*t519*t523*12-a1*mu*t523*t533*4+a3*mu*t519*t523*4-la*s2*s3*t523*5-mu*s2*s3*t523*4-la*t368*t519*t523*4-la*t518*t521*t533*3-mu*t518*t521*t533*4-a*a1*la*s3*t523*9+a*a3*la*s1*t523*9-a*a1*mu*s3*t523*3+a*a3*mu*s1*t523*3-a1*la*s1*s3*t523*12-a1*la*s2*s3*t523*12+a3*la*s1*s3*t523*12+a1*k*t519*t523*t533*4+a1*k*t521*t523*t533*4-a3*k*t518*t519*t523*4-a3*k*t519*t521*t523*4-a3*k*t519*t523*t533*4-a*la*s1*t368*t523*2+a*la*s1*t518*t519+a*la*s1*t518*t521*2-a*la*s1*t518*t533-a*la*s1*t521*t533*3+a*la*s2*t519*t521*8-a*la*s2*t521*t533*8-a*la*s3*t518*t521*2+a*la*s3*t519*t521*3-a*la*s3*t518*t533-a*la*s3*t521*t533*2-a1*mu*s1*s3*t523*4-a1*mu*s2*s3*t523*4+a3*mu*s1*s3*t523*4+a*mu*s1*t518*t521*2-a*mu*s1*t518*t533*2+a*mu*s2*t519*t521*8-a*mu*s2*t521*t533*8-a*mu*s3*t518*t521*2+a1*la*t368*t523*t533*4-a1*la*t518*t521*t533*4-a1*la*t519*t521*t533*4-a3*la*t368*t519*t523*4+a3*la*t519*t521*t533*4-a1*mu*t518*t521*t533*2-a1*mu*t519*t521*t533*2-a1*mu*t518*t523*t533*6-a1*mu*t519*t523*t533*6-a1*mu*t521*t523*t533*6+a3*mu*t519*t521*t523*6+a3*mu*t519*t521*t533*2+a3*mu*t519*t523*t533*6-la*s1*s2*t368*t523*4+la*s1*s2*t519*t521*3-la*s1*s2*t521*t533*4+la*s1*s3*t519*t521-la*s1*s3*t521*t533-la*s2*s3*t521*t533*3-mu*s1*s2*t521*t533*4+mu*s1*s3*t518*t519+mu*s1*s3*t519*t521-mu*s1*s3*t518*t533-mu*s1*s3*t521*t533-la*t368*t518*t519*t521*4-la*s1*s2*t368*t519*t521*4-la*s1*s3*t368*t518*t519-la*s1*s3*t368*t519*t521+la*s1*s3*t368*t518*t533+la*s1*s3*t368*t521*t533-la*s2*s3*t368*t519*t521*4+la*s2*s3*t368*t521*t533*4+a*a1*k*s1*t523*t533*8+a*a1*k*s2*t523*t533*8+a*a1*k*s3*t518*t523*2-a*a3*k*s1*t518*t523*2+a*a1*k*s3*t519*t523*2-a*a3*k*s1*t519*t523*2+a*a1*k*s3*t521*t523*2-a*a3*k*s1*t521*t523*2+a*a1*k*s3*t523*t533*2-a*a3*k*s1*t523*t533*2-a*a3*k*s2*t519*t523*8-a*a3*k*s3*t519*t523*8-a*a1*la*s1*t518*t533*2-a*a1*la*s1*t521*t533*14+a*a1*la*s3*t368*t523*3-a*a3*la*s1*t368*t523*3-a*a1*la*s2*t519*t533*2-a*a1*la*s2*t521*t533*14-a*a1*la*s3*t518*t519+a*a3*la*s1*t518*t519-a*a1*la*s3*t518*t521*4+a*a3*la*s1*t518*t521*4-a*a1*la*s3*t519*t521*4+a*a3*la*s1*t519*t521*4-a*a1*la*s3*t518*t533+a*a3*la*s1*t518*t533-a*a1*la*s3*t519*t533+a*a3*la*s1*t519*t533-a*a1*la*s3*t521*t533*4+a*a3*la*s1*t521*t533*4+a*a3*la*s2*t519*t521*14+a*a3*la*s2*t519*t533*2+a*a3*la*s3*t519*t521*14-a*a1*mu*s1*t521*t533*4-a*a1*mu*s1*t523*t533*12-a*a1*mu*s2*t521*t533*4-a*a1*mu*s2*t523*t533*12-a*a1*mu*s3*t518*t519+a*a3*mu*s1*t518*t519-a*a1*mu*s3*t518*t521*2+a*a3*mu*s1*t518*t521*2-a*a1*mu*s3*t519*t521*2+a*a3*mu*s1*t519*t521*2-a*a1*mu*s3*t518*t523*3+a*a3*mu*s1*t518*t523*3-a*a1*mu*s3*t519*t523*3+a*a3*mu*s1*t519*t523*3-a*a1*mu*s3*t521*t523*3+a*a3*mu*s1*t521*t523*3-a*a1*mu*s3*t518*t533+a*a3*mu*s1*t518*t533-a*a1*mu*s3*t519*t533+a*a3*mu*s1*t519*t533-a*a1*mu*s3*t521*t533*2+a*a3*mu*s1*t521*t533*2-a*a1*mu*s3*t523*t533*3+a*a3*mu*s1*t523*t533*3+a*a3*mu*s2*t519*t521*4+a*a3*mu*s2*t519*t523*12+a*a3*mu*s3*t519*t521*4+a*a3*mu*s3*t519*t523*12+a1*k*s1*s2*t523*t533*16+a1*k*s1*s3*t518*t523*4+a1*k*s1*s3*t521*t523*4+a1*k*s1*s3*t523*t533*4-a3*k*s1*s2*t519*t523*4-a3*k*s1*s2*t521*t523*4+a1*k*s2*s3*t523*t533*4-a3*k*s1*s2*t523*t533*4-a3*k*s1*s3*t518*t523*4-a3*k*s1*s3*t519*t523*4-a3*k*s1*s3*t521*t523*4-a3*k*s2*s3*t519*t523*16-a*la*s1*s2*s3*t533-a*mu*s1*s2*s3*t533*2+a1*la*s1*s3*t368*t523*4-a1*la*s1*s2*t521*t533*12-a3*la*s1*s2*t368*t523*4-a1*la*s1*s3*t518*t521*4-a1*la*s1*s3*t521*t533*4-a3*la*s1*s3*t368*t523*4-a1*la*s2*s3*t519*t521*4+a3*la*s1*s2*t519*t521*4-a1*la*s2*s3*t521*t533*4+a3*la*s1*s3*t518*t521*4+a3*la*s1*s3*t519*t521*4+a3*la*s2*s3*t519*t521*12-a3*k*t518*t519*t521*t533*4-a*la*s1*t368*t518*t519*2-a*la*s1*t368*t519*t521*2+a*la*s1*t368*t521*t533*2-a*la*s2*t368*t519*t521*8+a*la*s2*t368*t521*t533*8-a*la*s3*t368*t518*t519*2-a*la*s3*t368*t519*t521*2+a*la*s3*t368*t518*t533*2-a1*mu*s1*s2*t523*t533*24-a1*mu*s1*s3*t518*t521*2-a1*mu*s1*s3*t518*t523*6-a1*mu*s1*s3*t521*t523*6-a1*mu*s1*s3*t521*t533*2-a1*mu*s1*s3*t523*t533*6-a1*mu*s2*s3*t519*t521*2+a3*mu*s1*s2*t519*t521*2-a1*mu*s2*s3*t519*t523*6+a3*mu*s1*s2*t519*t523*6-a1*mu*s2*s3*t521*t523*6-a1*mu*s2*s3*t521*t533*2-a1*mu*s2*s3*t523*t533*6+a3*mu*s1*s3*t518*t521*2+a3*mu*s1*s3*t519*t521*2+a3*mu*s1*s3*t518*t523*6+a3*mu*s1*s3*t519*t523*6+a3*mu*s1*s3*t521*t523*6+a3*mu*s2*s3*t519*t523*24+a1*la*t368*t519*t521*t533*2-a3*la*t368*t518*t519*t521*2-a3*la*t368*t519*t521*t533*2-a1*mu*t518*t519*t521*t533*6-a*a3*k*s1*s2*s3*t523*8-a*a1*la*s1*s2*s3*t521*14-a*a1*la*s1*s2*s3*t533*2+a*a3*la*s1*s2*s3*t519*2+a*a1*k*s2*t519*t521*t533*8+a*a1*k*s3*t518*t519*t521*2-a*a3*k*s1*t518*t519*t521*2+a*a1*k*s3*t518*t519*t533*2-a*a3*k*s1*t518*t519*t533*2+a*a1*k*s3*t518*t521*t533*2-a*a3*k*s1*t518*t521*t533*2+a*a1*k*s3*t519*t521*t533*2-a*a3*k*s1*t519*t521*t533*2-a*a3*k*s2*t519*t521*t533*8-a*a3*k*s3*t518*t519*t521*8-a*a1*mu*s1*s2*s3*t521*4-a*a1*mu*s1*s2*s3*t523*12+a*a1*la*s1*t368*t521*t533*4+a*a1*la*s2*t368*t521*t533*4+a*a1*la*s3*t368*t518*t519-a*a3*la*s1*t368*t518*t519+a*a1*la*s3*t368*t518*t521*2-a*a3*la*s1*t368*t518*t521*2+a*a1*la*s3*t368*t519*t521*2-a*a3*la*s1*t368*t519*t521*2+a*a1*la*s3*t368*t518*t533-a*a3*la*s1*t368*t518*t533+a*a1*la*s3*t368*t519*t533-a*a3*la*s1*t368*t519*t533+a*a1*la*s3*t368*t521*t533*2-a*a3*la*s1*t368*t521*t533*2-a*a3*la*s2*t368*t519*t521*4-a*a3*la*s3*t368*t519*t521*4-a*a1*mu*s1*t518*t521*t533*12-a*a1*mu*s2*t519*t521*t533*12-a*a1*mu*s3*t518*t519*t521*3+a*a3*mu*s1*t518*t519*t521*3-a*a1*mu*s3*t518*t519*t533*3+a*a3*mu*s1*t518*t519*t533*3-a*a1*mu*s3*t518*t521*t533*3+a*a3*mu*s1*t518*t521*t533*3-a*a1*mu*s3*t519*t521*t533*3+a*a3*mu*s1*t519*t521*t533*3+a*a3*mu*s2*t519*t521*t533*12+a1*k*s1*s3*t518*t521*t533*4+a1*k*s2*s3*t519*t521*t533*4-a3*k*s1*s2*t519*t521*t533*4-a3*k*s1*s3*t518*t519*t521*4-a*la*s1*s2*s3*t368*t519*2+a1*la*s1*s3*t368*t518*t521*2+a1*la*s1*s3*t368*t521*t533*2-a3*la*s1*s2*t368*t519*t521*2+a1*la*s2*s3*t368*t521*t533*2-a3*la*s1*s2*t368*t521*t533*2-a3*la*s1*s3*t368*t518*t521*2-a3*la*s1*s3*t368*t519*t521*2-a1*mu*s1*s3*t518*t521*t533*6-a1*mu*s2*s3*t519*t521*t533*6+a3*mu*s1*s2*t519*t521*t533*6+a3*mu*s1*s3*t518*t519*t521*6+a*a1*k*s1*s2*s3*t521*t533*8-a*a3*k*s1*s2*s3*t519*t521*8-a*a3*la*s1*s2*s3*t368*t521*4-a*a1*mu*s1*s2*s3*t521*t533*12+a*a3*mu*s1*s2*s3*t519*t521*12))/(t519-t533);
    dP_dF.x3223 = (t364*t374*t386*(t584-t587-t592-t597+t618-t619+t620-t621-t622+t623+t624-t625+t626-t627-t628-t629+t630+t631+t632-t633-t634+t635-t636+t637-t638+t639-t640+t641-t642+t643-t644+t645-t646+t647-t648+t649-t650+t651-t652+t653-t654-t655+t656+t657+t658-t659-t660+t661-t662-t663+t664+t665-t666-t667+t668+t669+t670-t671+la*t518*t523*5-la*t523*t533*5+mu*t518*t523*4-mu*t523*t533*4+a*la*s2*t523*3-a*la*s3*t523*3+a*mu*s2*t523*2-a*mu*s3*t523*2-a2*la*t523*t533*12+a3*la*t518*t523*12-a2*mu*t523*t533*4+a3*mu*t518*t523*4-a*a2*la*s3*t523*9+a*a3*la*s2*t523*9-a*a2*mu*s3*t523*3+a*a3*mu*s2*t523*3-a2*la*s2*s3*t523*12+a3*la*s2*s3*t523*12+a2*k*t518*t523*t533*4+a2*k*t519*t523*t533*4+a2*k*t521*t523*t533*4-a3*k*t518*t519*t523*4-a3*k*t518*t521*t523*4-a3*k*t518*t523*t533*4+a*la*s1*t518*t521*8-a*la*s1*t521*t533*8+a*la*s2*t518*t519+a*la*s2*t518*t521*2+a*la*s2*t519*t521*2-a*la*s2*t521*t533*3+a*la*s3*t518*t521*3-a*la*s3*t519*t521*2-a*la*s3*t519*t533-a*la*s3*t521*t533*2-a2*mu*s2*s3*t523*4+a3*mu*s2*s3*t523*4+a*mu*s1*t518*t521*8-a*mu*s1*t521*t533*8+a*mu*s2*t518*t519*2+a*mu*s2*t518*t521*2+a*mu*s3*t518*t521*2-a*mu*s3*t519*t533*2+a2*la*t368*t523*t533*4-a3*la*t368*t518*t523*4-a2*la*t518*t521*t533*4+a3*la*t518*t521*t533*4-a2*mu*t518*t521*t533*2-a2*mu*t518*t523*t533*6-a2*mu*t521*t523*t533*6+a3*mu*t518*t521*t523*6+a3*mu*t518*t521*t533*2+a3*mu*t518*t523*t533*6-la*s1*s2*t368*t523*4+la*s1*s3*t368*t523*4+la*s1*s2*t518*t521*3-la*s1*s2*t521*t533*4+la*s1*s3*t518*t521*4-la*s1*s3*t521*t533*3+la*s2*s3*t518*t521-la*s2*s3*t521*t533+mu*s1*s2*t518*t521*4-mu*s1*s2*t521*t533*4+mu*s1*s3*t518*t521*4-mu*s1*s3*t521*t533*4+mu*s2*s3*t518*t519+mu*s2*s3*t518*t521-mu*s2*s3*t519*t533-mu*s2*s3*t521*t533-la*t368*t518*t519*t521*4+la*t368*t519*t521*t533*4-la*s1*s2*t368*t518*t521*4+la*s1*s3*t368*t521*t533*4-la*s2*s3*t368*t518*t519-la*s2*s3*t368*t518*t521+la*s2*s3*t368*t519*t533+la*s2*s3*t368*t521*t533+a*a2*k*s1*t523*t533*8-a*a3*k*s1*t518*t523*8+a*a2*k*s2*t523*t533*8+a*a2*k*s3*t518*t523*2-a*a3*k*s2*t518*t523*2+a*a2*k*s3*t519*t523*2-a*a3*k*s2*t519*t523*2+a*a2*k*s3*t521*t523*2-a*a3*k*s2*t521*t523*2+a*a2*k*s3*t523*t533*2-a*a3*k*s2*t523*t533*2-a*a3*k*s3*t518*t523*8-a*a2*la*s1*t518*t533*2-a*a2*la*s1*t521*t533*14+a*a2*la*s3*t368*t523*3-a*a3*la*s2*t368*t523*3+a*a3*la*s1*t518*t521*14+a*a3*la*s1*t518*t533*2-a*a2*la*s2*t521*t533*14-a*a2*la*s3*t518*t519+a*a3*la*s2*t518*t519-a*a2*la*s3*t518*t521*4+a*a3*la*s2*t518*t521*4-a*a2*la*s3*t519*t521*4+a*a3*la*s2*t519*t521*4-a*a2*la*s3*t518*t533+a*a3*la*s2*t518*t533-a*a2*la*s3*t519*t533+a*a3*la*s2*t519*t533-a*a2*la*s3*t521*t533*4+a*a3*la*s2*t521*t533*4+a*a3*la*s3*t518*t521*14-a*a2*mu*s1*t521*t533*4-a*a2*mu*s1*t523*t533*12+a*a3*mu*s1*t518*t521*4+a*a3*mu*s1*t518*t523*12-a*a2*mu*s2*t521*t533*4-a*a2*mu*s2*t523*t533*12-a*a2*mu*s3*t518*t519+a*a3*mu*s2*t518*t519-a*a2*mu*s3*t518*t521*2+a*a3*mu*s2*t518*t521*2-a*a2*mu*s3*t519*t521*2+a*a3*mu*s2*t519*t521*2-a*a2*mu*s3*t518*t523*3+a*a3*mu*s2*t518*t523*3-a*a2*mu*s3*t519*t523*3+a*a3*mu*s2*t519*t523*3-a*a2*mu*s3*t521*t523*3+a*a3*mu*s2*t521*t523*3-a*a2*mu*s3*t518*t533+a*a3*mu*s2*t518*t533-a*a2*mu*s3*t519*t533+a*a3*mu*s2*t519*t533-a*a2*mu*s3*t521*t533*2+a*a3*mu*s2*t521*t533*2-a*a2*mu*s3*t523*t533*3+a*a3*mu*s2*t523*t533*3+a*a3*mu*s3*t518*t521*4+a*a3*mu*s3*t518*t523*12+a2*k*s1*s2*t523*t533*16+a2*k*s1*s3*t518*t523*4-a3*k*s1*s2*t518*t523*4+a2*k*s1*s3*t521*t523*4-a3*k*s1*s2*t521*t523*4+a2*k*s1*s3*t523*t533*4-a3*k*s1*s2*t523*t533*4-a3*k*s1*s3*t518*t523*16+a2*k*s2*s3*t519*t523*4+a2*k*s2*s3*t521*t523*4+a2*k*s2*s3*t523*t533*4-a3*k*s2*s3*t518*t523*4-a3*k*s2*s3*t519*t523*4-a3*k*s2*s3*t521*t523*4+a*la*s1*s2*s3*t518-a*la*s1*s2*s3*t533+a*mu*s1*s2*s3*t518*2-a*mu*s1*s2*s3*t533*2+a2*la*s1*s3*t368*t523*4-a3*la*s1*s2*t368*t523*4-a2*la*s1*s2*t521*t533*12+a2*la*s2*s3*t368*t523*4+a3*la*s1*s2*t518*t521*4-a2*la*s1*s3*t521*t533*4-a3*la*s2*s3*t368*t523*4+a3*la*s1*s3*t518*t521*12-a2*la*s2*s3*t519*t521*4-a2*la*s2*s3*t521*t533*4+a3*la*s2*s3*t518*t521*4+a3*la*s2*s3*t519*t521*4+a2*k*t518*t519*t521*t533*4-a3*k*t518*t519*t521*t533*4-a*la*s1*t368*t518*t521*8+a*la*s1*t368*t521*t533*8-a*la*s2*t368*t518*t519*2-a*la*s2*t368*t519*t521*2+a*la*s2*t368*t519*t533*2-a*la*s3*t368*t518*t519*2+a*la*s3*t368*t519*t521*2+a*la*s3*t368*t519*t533*2-a2*mu*s1*s2*t523*t533*24+a3*mu*s1*s2*t518*t521*2+a3*mu*s1*s2*t518*t523*6-a2*mu*s1*s3*t521*t533*2-a2*mu*s1*s3*t523*t533*6-a2*mu*s2*s3*t519*t521*2+a3*mu*s1*s3*t518*t523*24-a2*mu*s2*s3*t519*t523*6-a2*mu*s2*s3*t521*t523*6-a2*mu*s2*s3*t521*t533*2-a2*mu*s2*s3*t523*t533*6+a3*mu*s2*s3*t518*t521*2+a3*mu*s2*s3*t519*t521*2+a3*mu*s2*s3*t518*t523*6+a3*mu*s2*s3*t519*t523*6+a3*mu*s2*s3*t521*t523*6+a2*la*t368*t518*t521*t533*2+a2*la*t368*t519*t521*t533*2-a3*la*t368*t518*t519*t521*2-a3*la*t368*t518*t521*t533*2+a*a2*k*s1*s2*s3*t523*8-a*a3*k*s1*s2*s3*t523*8-a*a2*la*s1*s2*s3*t533*2+a*a3*la*s1*s2*s3*t518*2+a*a2*k*s1*t518*t521*t533*8-a*a3*k*s1*t518*t521*t533*8+a*a2*k*s2*t519*t521*t533*8+a*a2*k*s3*t518*t519*t521*2-a*a3*k*s2*t518*t519*t521*2+a*a2*k*s3*t518*t519*t533*2-a*a3*k*s2*t518*t519*t533*2+a*a2*k*s3*t518*t521*t533*2-a*a3*k*s2*t518*t521*t533*2+a*a2*k*s3*t519*t521*t533*2-a*a3*k*s2*t519*t521*t533*2-a*a3*k*s3*t518*t519*t521*8+a*a2*la*s1*t368*t521*t533*4-a*a3*la*s1*t368*t518*t521*4+a*a2*la*s2*t368*t521*t533*4+a*a2*la*s3*t368*t518*t519-a*a3*la*s2*t368*t518*t519+a*a2*la*s3*t368*t518*t521*2-a*a3*la*s2*t368*t518*t521*2+a*a2*la*s3*t368*t519*t521*2-a*a3*la*s2*t368*t519*t521*2+a*a2*la*s3*t368*t518*t533-a*a3*la*s2*t368*t518*t533+a*a2*la*s3*t368*t519*t533-a*a3*la*s2*t368*t519*t533+a*a2*la*s3*t368*t521*t533*2-a*a3*la*s2*t368*t521*t533*2-a*a3*la*s3*t368*t518*t521*4-a*a2*mu*s1*t518*t521*t533*12+a*a3*mu*s1*t518*t521*t533*12-a*a2*mu*s3*t518*t519*t521*3+a*a3*mu*s2*t518*t519*t521*3-a*a2*mu*s3*t518*t519*t533*3+a*a3*mu*s2*t518*t519*t533*3-a*a2*mu*s3*t518*t521*t533*3+a*a3*mu*s2*t518*t521*t533*3-a*a2*mu*s3*t519*t521*t533*3+a*a3*mu*s2*t519*t521*t533*3+a2*k*s1*s3*t518*t521*t533*4-a3*k*s1*s2*t518*t521*t533*4+a2*k*s2*s3*t519*t521*t533*4-a3*k*s2*s3*t518*t519*t521*4+a2*la*s1*s3*t368*t518*t521*2-a3*la*s1*s2*t368*t518*t521*2+a2*la*s1*s3*t368*t521*t533*2-a3*la*s1*s2*t368*t521*t533*2+a2*la*s2*s3*t368*t519*t521*2+a2*la*s2*s3*t368*t521*t533*2-a3*la*s2*s3*t368*t518*t521*2-a3*la*s2*s3*t368*t519*t521*2-a2*mu*s1*s3*t518*t521*t533*6+a3*mu*s1*s2*t518*t521*t533*6-a2*mu*s2*s3*t519*t521*t533*6+a3*mu*s2*s3*t518*t519*t521*6+a*a2*k*s1*s2*s3*t521*t533*8-a*a3*k*s1*s2*s3*t518*t521*8+a*a2*la*s1*s2*s3*t368*t521*4-a*a3*la*s1*s2*s3*t368*t521*4-a*a2*mu*s1*s2*s3*t521*t533*12+a*a3*mu*s1*s2*s3*t518*t521*12))/(t518-t533);
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
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
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
template class NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<float,2>;
template class NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<double,2>;
template class NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<double,3>;
#endif
