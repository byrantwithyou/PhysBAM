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
    //s1o=-1.0; s2o=-1.0; s3o=-1.0;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
~NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA()
{
}

//#####################################################################
// Get a by applying Newton's method and its partial derivatives
// through implicit differentiation
//#####################################################################
/*template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,3>::
 Calculate_A(const T s1,const T s2,const T s3)
 {
 /////////////Hopefully I can get the member function to work at some point
 int maxiter = 100; int iter;
 T tol = (T)1e-8;
 //if (d==2){PHYSBAM_FATAL_ERROR();}
 //if (pow(s1-s1o,2.0)+pow(s1-s1o,2.0)+pow(s1-s1o,2.0)<tol) {return;} //Don't need to recalculate
 //T c = extrapolation_cutoff;    
 T c13 = pow(c,(double)1/(double)3);
 T a0,f,fa;
 
 if (s1<0) {a0=c13-s1;}          //Pick initial guess
 else{ if (s2<0) {a0=c13-s2;}
 else{ if (s3<0) {a0=c13-s3;}
 else a0=c13;}}
 
 iter=0; f=2e-8; a=a0;
 
 while(f>tol){                   //Newton iteration
 f=(s1+a)*(s2+a)*(s3+a)-c;
 fa=(s1+a)*(s2+a)+(s1+a)*(s3+a)+(s3+a)*(s2+a);
 if (abs(fa)<tol){PHYSBAM_FATAL_ERROR();}
 a=a-f/fa;
 iter++; if (iter>maxiter){PHYSBAM_FATAL_ERROR();}
 }
 
 //Now we get a1,a2,a3,a11,a12,a13,a22,a23,a33
 
 T den = (3.0*a*a-2.0*a*(s1+s2+s3)+(s1*s2+s1*s3+s2*s3));
 
 a1 = -(a+s2)*(a+s3)/den;
 a2 = -(a+s1)*(a+s3)/den;
 a3 = -(a+s1)*(a+s2)/den;
 a11 = -2.0*a1*(s3+2.0*a+a1*s3+3.0*a1*a+s2+a1*s2+a1*s1)/den;
 a22 = -2.0*a2*(s3+2.0*a+a2*s3+3.0*a2*a+s1+a2*s1+a2*s2)/den;
 a33 = -2.0*a3*(s1+2.0*a+a3*s1+3.0*a3*a+s2+a3*s2+a3*s3)/den;
 a12 = -(s3+a+a1*s3+2.0*a1*a+a1*s1+2.0*a1*a2*s3+6.0*a1*a2*a+2.0*a1*a2*s2+2.0*a1*a2*s1+a2*s3+2.0*a2*a+a2*s2)/den;
 a13 = -(s2+a+a1*s2+2.0*a1*a+a1*s1+2.0*a1*a3*s2+6.0*a1*a3*a+2.0*a1*a3*s3+2.0*a1*a3*s1+a3*s2+2.0*a3*a+a3*s3)/den;
 a23 = -(s1+a+a3*s1+2.0*a3*a+a3*s3+2.0*a3*a2*s1+6.0*a3*a2*a+2.0*a3*a2*s2+2.0*a3*a2*s3+a2*s1+2.0*a2*a+a2*s2)/den;
 
 /////////////Hopefully I can get the member function to work at some point
 
 } */
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
    t3 = s2*(1.0/2.0);
    t4 = c*4.0;
    t5 = t2*t2;
    t6 = t4+t5;
    t7 = sqrt(t6);
    t8 = t7*(1.0/2.0);
    t10 = s1*(1.0/2.0);
    t9 = -t10+t3+t8;
    t11 = 1.0/t9;
    t12 = t10-t3+t8;
    t13 = t12*t9;
    t14 = log(t13);
    t15 = 1.0/t12;
    t16 = t10+t3-t8;
    t0 = -mu*t14+t16*(mu*t12-mu*t15+la*t14*t15)+t16*(-mu*t11+mu*t9+la*t11*t14)+la*(t14*t14)*(1.0/2.0)+mu*(t12*t12+t9*t9-2.0)*(1.0/2.0);
    
    
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
    //if (pow(s1-s1o,2.0)+pow(s1-s1o,2.0)+pow(s1-s1o,2.0)<tol) {return;} //Don't need to recalculate
    T c = extrapolation_cutoff;    
    T c13 = pow(c,(double)1/(double)3);
    T a0,f,fa;
    
    if (s1<0) {a0=c13-s1;}          //Pick initial guess
    else{ if (s2<0) {a0=c13-s2;}
    else{ if (s3<0) {a0=c13-s3;}
        else a0=c13;}}
    
    iter=0; f=2e-8; a=a0;
    
    while(f>tol){                   //Newton iteration
        f=(s1+a)*(s2+a)*(s3+a)-c;
        fa=(s1+a)*(s2+a)+(s1+a)*(s3+a)+(s3+a)*(s2+a);
        if (abs(fa)<tol){PHYSBAM_FATAL_ERROR();}
        a=a-f/fa;
        iter++; if (iter>maxiter){PHYSBAM_FATAL_ERROR();}
    }
    
    //Now we get a1,a2,a3,a11,a12,a13,a22,a23,a33
    
    T den = (3.0*a*a-2.0*a*(s1+s2+s3)+(s1*s2+s1*s3+s2*s3));
    
    a1 = -(a+s2)*(a+s3)/den;
    a2 = -(a+s1)*(a+s3)/den;
    a3 = -(a+s1)*(a+s2)/den;
    a11 = -2.0*a1*(s3+2.0*a+a1*s3+3.0*a1*a+s2+a1*s2+a1*s1)/den;
    a22 = -2.0*a2*(s3+2.0*a+a2*s3+3.0*a2*a+s1+a2*s1+a2*s2)/den;
    a33 = -2.0*a3*(s1+2.0*a+a3*s1+3.0*a3*a+s2+a3*s2+a3*s3)/den;
    a12 = -(s3+a+a1*s3+2.0*a1*a+a1*s1+2.0*a1*a2*s3+6.0*a1*a2*a+2.0*a1*a2*s2+2.0*a1*a2*s1+a2*s3+2.0*a2*a+a2*s2)/den;
    a13 = -(s2+a+a1*s2+2.0*a1*a+a1*s1+2.0*a1*a3*s2+6.0*a1*a3*a+2.0*a1*a3*s3+2.0*a1*a3*s1+a3*s2+2.0*a3*a+a3*s3)/den;
    a23 = -(s1+a+a3*s1+2.0*a3*a+a3*s3+2.0*a3*a2*s1+6.0*a3*a2*a+2.0*a3*a2*s2+2.0*a3*a2*s3+a2*s1+2.0*a2*a+a2*s2)/den;
    
    /////////////Hopefully I can get the member function to work at some point   
    
    T mu = constant_mu;
    T la = constant_lambda;
    //Calculate_A(s1,s2,s3);
    T k = extra_force_coefficient; 
    
    T t0 = -mu*log((a+s1)*(a+s2)*(a+s3))-a*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+(a*a)*k+la*pow(log((a+s1)*(a+s2)*(a+s3)),2.0)*(1.0/2.0)+mu*(pow(a+s1,2.0)+pow(a+s2,2.0)+pow(a+s3,2.0)-2.0)*(1.0/2.0);
    
    
    
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
    t19 = c*4.0;
    t20 = t18*t18;
    t21 = t19+t20;
    t22 = s1*2.0;
    t23 = s2*2.0;
    t24 = t22-t23;
    t25 = 1.0/sqrt(t21);
    t26 = t24*t25*(1.0/4.0);
    t27 = t26-1.0/2.0;
    t28 = s1*(1.0/2.0);
    t29 = s2*(1.0/2.0);
    t30 = sqrt(t21);
    t31 = t30*(1.0/2.0);
    t32 = -t28+t29+t31;
    t33 = 1.0/(t32*t32);
    t34 = t28-t29+t31;
    t35 = t26+1.0/2.0;
    t36 = t27*t34;
    t37 = t32*t35;
    t38 = t36+t37;
    t39 = 1.0/(t34*t34);
    t40 = t32*t34;
    t41 = log(t40);
    t42 = t24*t25*(1.0/2.0);
    t43 = 1.0/t32;
    t44 = 1.0/t34;
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
    t55 = t42-1.0;
    t56 = t32*t55;
    t57 = t42+1.0;
    t58 = t34*t57;
    t59 = t56+t58;
    t60 = mu*t59*(1.0/2.0);
    t61 = mu*t32;
    t62 = la*t41*t43;
    t63 = t61+t62-mu*t43;
    t64 = mu*t34;
    t65 = la*t41*t44;
    t66 = t64+t65-mu*t44;
    t67 = la*t38*t41*t43*t44;
    result.x11 = t54+t60+t67-t27*t63-t27*t66+t49*(t28+t29-t30*(1.0/2.0))-mu*t38*t43*t44;
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
    //if (d==2){PHYSBAM_FATAL_ERROR();}
    //if (pow(s1-s1o,2.0)+pow(s1-s1o,2.0)+pow(s1-s1o,2.0)<tol) {return;} //Don't need to recalculate
    T c = extrapolation_cutoff;    
    T c13 = pow(c,(double)1/(double)3);
    T a0,f,fa;
    
    if (s1<0) {a0=c13-s1;}          //Pick initial guess
    else{ if (s2<0) {a0=c13-s2;}
    else{ if (s3<0) {a0=c13-s3;}
        else a0=c13;}}
    
    iter=0; f=2e-8; a=a0;
    
    while(f>tol){                   //Newton iteration
        f=(s1+a)*(s2+a)*(s3+a)-c;
        fa=(s1+a)*(s2+a)+(s1+a)*(s3+a)+(s3+a)*(s2+a);
        if (abs(fa)<tol){PHYSBAM_FATAL_ERROR();}
        a=a-f/fa;
        iter++; if (iter>maxiter){PHYSBAM_FATAL_ERROR();}
    }
    
    //Now we get a1,a2,a3,a11,a12,a13,a22,a23,a33
    
    T den = (3.0*a*a-2.0*a*(s1+s2+s3)+(s1*s2+s1*s3+s2*s3));
    
    a1 = -(a+s2)*(a+s3)/den;
    a2 = -(a+s1)*(a+s3)/den;
    a3 = -(a+s1)*(a+s2)/den;
    a11 = -2.0*a1*(s3+2.0*a+a1*s3+3.0*a1*a+s2+a1*s2+a1*s1)/den;
    a22 = -2.0*a2*(s3+2.0*a+a2*s3+3.0*a2*a+s1+a2*s1+a2*s2)/den;
    a33 = -2.0*a3*(s1+2.0*a+a3*s1+3.0*a3*a+s2+a3*s2+a3*s3)/den;
    a12 = -(s3+a+a1*s3+2.0*a1*a+a1*s1+2.0*a1*a2*s3+6.0*a1*a2*a+2.0*a1*a2*s2+2.0*a1*a2*s1+a2*s3+2.0*a2*a+a2*s2)/den;
    a13 = -(s2+a+a1*s2+2.0*a1*a+a1*s1+2.0*a1*a3*s2+6.0*a1*a3*a+2.0*a1*a3*s3+2.0*a1*a3*s1+a3*s2+2.0*a3*a+a3*s3)/den;
    a23 = -(s1+a+a3*s1+2.0*a3*a+a3*s3+2.0*a3*a2*s1+6.0*a3*a2*a+2.0*a3*a2*s2+2.0*a3*a2*s3+a2*s1+2.0*a2*a+a2*s2)/den;
    
    /////////////Hopefully I can get the member function to work at some point
    
    T mu = constant_mu;
    T la = constant_lambda;
    DIAGONAL_MATRIX<T,3> result;
    //Calculate_A(s1,s2,s3);
    T k = extra_force_coefficient;    
    /*    T t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42;
     
     
     t2 = a1+1.0;
     t3 = a+s1;
     t4 = 1.0/(t3*t3);
     t5 = a+s2;
     t6 = a+s3;
     t7 = t3*t5*t6;
     t8 = log(t7);
     t9 = 1.0/t3;
     t10 = 1.0/t5;
     t11 = a*2.0;
     t12 = 1.0/(t5*t5);
     t13 = 1.0/t6;
     t14 = a1*t3*t5;
     t15 = a1*t3*t6;
     t16 = t2*t5*t6;
     t17 = t14+t15+t16;
     t18 = a2+1.0;
     t19 = mu*t3;
     t20 = la*t8*t9;
     t35 = mu*t9;
     t21 = t19+t20-t35;
     t22 = mu*t5;
     t23 = la*t10*t8;
     t36 = mu*t10;
     t24 = t22+t23-t36;
     t25 = s2*2.0;
     t26 = t11+t25;
     t27 = s1*2.0;
     t28 = t11+t27;
     t29 = s3*2.0;
     t30 = t11+t29;
     t31 = a2*t3*t5;
     t32 = a2*t5*t6;
     t33 = t18*t3*t6;
     t34 = t31+t32+t33;
     t37 = a3+1.0;
     t38 = a3*mu;
     t39 = a3*t3*t6;
     t40 = a3*t5*t6;
     t41 = t3*t37*t5;
     t42 = t39+t40+t41;
     result.x11 = -a1*t21-a1*t24-a*(a1*mu+a1*mu*t12-a1*la*t12*t8+la*t12*t13*t17*t9)-a*(mu*t2+mu*t2*t4-la*t2*t4*t8+la*t10*t13*t17*t4)+mu*(a1*t26+a1*t30+t2*t28)*(1.0/2.0)-mu*t10*t13*t17*t9+la*t10*t13*t17*t8*t9+2.0*k*a1*a;
     result.x22 = -a2*t21-a2*t24-a*(a2*mu+a2*mu*t4-a2*la*t4*t8+la*t10*t13*t34*t4)-a*(mu*t18+mu*t12*t18-la*t12*t18*t8+la*t12*t13*t34*t9)+mu*(a2*t28+a2*t30+t18*t26)*(1.0/2.0)-mu*t10*t13*t34*t9+la*t10*t13*t34*t8*t9+2.0*k*a2*a;
     result.x33 = -a3*t21-a3*t24+mu*(a3*t26+a3*t28+t30*t37)*(1.0/2.0)-a*(t38+a3*mu*t12-a3*la*t12*t8+la*t12*t13*t42*t9)-a*(t38+a3*mu*t4-a3*la*t4*t8+la*t10*t13*t4*t42)-mu*t10*t13*t42*t9+la*t10*t13*t42*t8*t9+2.0*k*a3*a;*/
    result.x11 = -a*(mu*(a1+1.0)+mu*1.0/pow(a+s1,2.0)*(a1+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a1+1.0)+(la*1.0/pow(a+s1,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3)))-a1*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a1*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a1*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+mu*((a1+1.0)*(a*2.0+s1*2.0)+a1*(a*2.0+s2*2.0)+a1*(a*2.0+s3*2.0))*(1.0/2.0)-a*(a1*mu+a1*mu*1.0/pow(a+s2,2.0)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+(la*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3)))-a*(a1*mu+a1*mu*1.0/pow(a+s3,2.0)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+(la*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2)))+a*a1*k*2.0-(mu*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2)*(a+s3));
    result.x22 = -a*(mu*(a2+1.0)+mu*1.0/pow(a+s2,2.0)*(a2+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a2+1.0)+(la*1.0/pow(a+s2,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3)))-a2*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a2*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a2*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+mu*((a2+1.0)*(a*2.0+s2*2.0)+a2*(a*2.0+s1*2.0)+a2*(a*2.0+s3*2.0))*(1.0/2.0)-a*(a2*mu+a2*mu*1.0/pow(a+s1,2.0)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)+(la*1.0/pow(a+s1,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3)))-a*(a2*mu+a2*mu*1.0/pow(a+s3,2.0)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+(la*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2)))+a*a2*k*2.0-(mu*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2)*(a+s3));
    result.x33 = -a*(mu*(a3+1.0)+mu*1.0/pow(a+s3,2.0)*(a3+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a3+1.0)+(la*1.0/pow(a+s3,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2)))-a3*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a3*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a3*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+mu*((a3+1.0)*(a*2.0+s3*2.0)+a3*(a*2.0+s1*2.0)+a3*(a*2.0+s2*2.0))*(1.0/2.0)-a*(a3*mu+a3*mu*1.0/pow(a+s1,2.0)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)+(la*1.0/pow(a+s1,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s2)*(a+s3)))-a*(a3*mu+a3*mu*1.0/pow(a+s2,2.0)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+(la*1.0/pow(a+s2,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s3)))+a*a3*k*2.0-(mu*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2)*(a+s3));
    
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
    t70 = c*4.0;
    t71 = t69*t69;
    t72 = t70+t71;
    t73 = pow(t72,3.0/2.0);
    t74 = s1*s1;
    t75 = s2*s2;
    t76 = t74*t74;
    t77 = t75*t75;
    t78 = sqrt(t72);
    t79 = s2*(1.0/2.0);
    t80 = t78*(1.0/2.0);
    t81 = s1*(1.0/2.0);
    t82 = -t79+t80+t81;
    t83 = t79+t80-t81;
    t84 = t82*t83;
    t85 = log(t84);
    t86 = t74*t75*1.2E1;
    t87 = s1*t73*2.0;
    t88 = c*t74*1.6E1;
    t89 = c*t75*1.6E1;
    t90 = c*c;
    t91 = t90*3.2E1;
    t92 = t76*2.0;
    t93 = t77*2.0;
    t94 = c*mu*3.2E1;
    t95 = mu*t76*2.0;
    t96 = mu*t77*2.0;
    t97 = mu*t74*t75*1.2E1;
    t98 = la*s1*s2*t85*1.6E1;
    t99 = la*s1*t78*t85*4.0;
    t100 = la*s2*t78*t85*4.0;
    t101 = mu*s1*t75*t78;
    t102 = mu*s2*t74*t78;
    t104 = s1*s2*2.0;
    t103 = -t104+t70+t74+t75;
    t105 = pow(t103,3.0/2.0);
    t106 = log(c);
    t107 = 1.0/c;
    t108 = sqrt(t103);
    t109 = s1+s2;if(fabs(t109)<panic_threshold) t109=t109<0?-panic_threshold:panic_threshold;
    t110 = 1.0/t109;
    t111 = 1.0/sqrt(t103);
    t112 = mu*t90*2.0;
    t113 = mu*s1*t108;
    t114 = mu*s2*t108;
    t115 = c*la*t106*2.0;
    dP_dF.x1111 = (t100+t101+t102+t94+t95+t96+t97+t98+t99+mu*t74*4.0+mu*t75*1.2E1-la*t74*log(t82*(s1*(-1.0/2.0)+t79+t80))*4.0-c*la*t85*3.2E1+c*mu*t74*1.2E1+c*mu*t75*4.0-mu*s1*s2*1.6E1-la*t75*t85*1.2E1+mu*s1*t73*3.0-mu*s1*t78*4.0-mu*s2*t73-mu*s2*t78*4.0-c*mu*s1*s2*1.6E1-mu*s1*s2*t74*8.0-mu*s1*s2*t75*8.0-mu*s1*t74*t78-mu*s2*t75*t78)/(t86+t87+t88+t89+t91+t92+t93-s2*t73*2.0-c*s1*s2*3.2E1-s1*s2*t74*8.0-s1*s2*t75*8.0);
    dP_dF.x2222 = (t100+t101+t102+t94+t95+t96+t97+t98+t99+mu*t74*1.2E1+mu*t75*4.0-c*la*t85*3.2E1+c*mu*t74*4.0+c*mu*t75*1.2E1-mu*s1*s2*1.6E1-la*t74*t85*1.2E1-la*t75*t85*4.0-mu*s1*t73-mu*s1*t78*4.0+mu*s2*t73*3.0-mu*s2*t78*4.0-c*mu*s1*s2*1.6E1-mu*s1*s2*t74*8.0-mu*s1*s2*t75*8.0-mu*s1*t74*t78-mu*s2*t75*t78)/(t86-t87+t88+t89+t91+t92+t93+s2*t73*2.0-c*s1*s2*3.2E1-s1*s2*t74*8.0-s1*s2*t75*8.0);
    dP_dF.x2211 = -1.0/pow(t103,3.0/2.0)*t107*(mu*t105-c*mu*s1*2.0-c*mu*s2*2.0-la*t105*t106+mu*s1*t90*2.0+mu*s2*t90*2.0+c*la*s1*t106*2.0+c*la*s2*t106*2.0);
    dP_dF.x2121 = t107*t110*t111*(t112+t113+t114+t115-c*mu*2.0-mu*t74-mu*t75+c*mu*t74+c*mu*t75+la*t106*t74+la*t106*t75-la*s1*t106*t108-la*s2*t106*t108);
    dP_dF.x2112 = -t107*t110*t111*(t112-t113-t114+t115-c*mu*2.0+mu*s1*s2*2.0-c*mu*s1*s2*2.0-la*s1*s2*t106*2.0+la*s1*t106*t108+la*s2*t106*t108);
    
    
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
    T k = extra_force_coefficient;
    
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
    /////////////Hopefully I can get the member function to work at some point
    T a,a1,a2,a3,a11,a12,a13,a22,a23,a33;
    int maxiter = 100; int iter;
    T tol = (T)1e-8;
    //if (d==2){PHYSBAM_FATAL_ERROR();}
    //if (pow(s1-s1o,2.0)+pow(s1-s1o,2.0)+pow(s1-s1o,2.0)<tol) {return;} //Don't need to recalculate
    //T c = extrapolation_cutoff;    
    T c13 = pow(c,(double)1/(double)3);
    T a0,f,fa;
    
    if (s1<0) {a0=c13-s1;}          //Pick initial guess
    else{ if (s2<0) {a0=c13-s2;}
    else{ if (s3<0) {a0=c13-s3;}
        else a0=c13;}}
    
    iter=0; f=2e-8; a=a0;
    
    while(f>tol){                   //Newton iteration
        f=(s1+a)*(s2+a)*(s3+a)-c;
        fa=(s1+a)*(s2+a)+(s1+a)*(s3+a)+(s3+a)*(s2+a);
        if (abs(fa)<tol){PHYSBAM_FATAL_ERROR();}
        a=a-f/fa;
        iter++; if (iter>maxiter){PHYSBAM_FATAL_ERROR();}
    }
    
    //Now we get a1,a2,a3,a11,a12,a13,a22,a23,a33
    
    T den = (3.0*a*a-2.0*a*(s1+s2+s3)+(s1*s2+s1*s3+s2*s3));
    
    a1 = -(a+s2)*(a+s3)/den;
    a2 = -(a+s1)*(a+s3)/den;
    a3 = -(a+s1)*(a+s2)/den;
    a11 = -2.0*a1*(s3+2.0*a+a1*s3+3.0*a1*a+s2+a1*s2+a1*s1)/den;
    a22 = -2.0*a2*(s3+2.0*a+a2*s3+3.0*a2*a+s1+a2*s1+a2*s2)/den;
    a33 = -2.0*a3*(s1+2.0*a+a3*s1+3.0*a3*a+s2+a3*s2+a3*s3)/den;
    a12 = -(s3+a+a1*s3+2.0*a1*a+a1*s1+2.0*a1*a2*s3+6.0*a1*a2*a+2.0*a1*a2*s2+2.0*a1*a2*s1+a2*s3+2.0*a2*a+a2*s2)/den;
    a13 = -(s2+a+a1*s2+2.0*a1*a+a1*s1+2.0*a1*a3*s2+6.0*a1*a3*a+2.0*a1*a3*s3+2.0*a1*a3*s1+a3*s2+2.0*a3*a+a3*s3)/den;
    a23 = -(s1+a+a3*s1+2.0*a3*a+a3*s3+2.0*a3*a2*s1+6.0*a3*a2*a+2.0*a3*a2*s2+2.0*a3*a2*s3+a2*s1+2.0*a2*a+a2*s2)/den;
    
    /////////////Hopefully I can get the member function to work at some point
    
    dP_dF.x1111 = a1*(mu*(a1+1.0)+mu*1.0/pow(a+s1,2.0)*(a1+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a1+1.0)+(la*1.0/pow(a+s1,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3)))*-2.0-a11*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a11*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a11*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+(a1*a1)*k*2.0+mu*(pow(a1+1.0,2.0)*2.0+a11*(a*2.0+s1*2.0)+a11*(a*2.0+s2*2.0)+a11*(a*2.0+s3*2.0)+(a1*a1)*4.0)*(1.0/2.0)+a*(-a11*mu+(a1*a1)*mu*1.0/pow(a+s2,3.0)*2.0-a11*mu*1.0/pow(a+s2,2.0)+a11*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)-(a1*a1)*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,3.0)*2.0-(la*1.0/pow(a+s2,2.0)*((a1*a1)*(a*2.0+s1*2.0)+a11*(a+s1)*(a+s2)+a11*(a+s1)*(a+s3)+a11*(a+s2)*(a+s3)+a1*(a1+1.0)*(a*2.0+s2*2.0)+a1*(a1+1.0)*(a*2.0+s3*2.0)))/((a+s1)*(a+s3))+(a1*la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s1)+(a1*la*1.0/pow(a+s2,3.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*3.0)/((a+s1)*(a+s3))+(la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a1+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s3))+a*(-a11*mu+(a1*a1)*mu*1.0/pow(a+s3,3.0)*2.0-a11*mu*1.0/pow(a+s3,2.0)+a11*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)-(a1*a1)*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,3.0)*2.0-(la*1.0/pow(a+s3,2.0)*((a1*a1)*(a*2.0+s1*2.0)+a11*(a+s1)*(a+s2)+a11*(a+s1)*(a+s3)+a11*(a+s2)*(a+s3)+a1*(a1+1.0)*(a*2.0+s2*2.0)+a1*(a1+1.0)*(a*2.0+s3*2.0)))/((a+s1)*(a+s2))+(a1*la*1.0/pow(a+s3,3.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*3.0)/((a+s1)*(a+s2))+(a1*la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s1)+(la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a1+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s2))-a1*(a1*mu+a1*mu*1.0/pow(a+s2,2.0)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+(la*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3)))*2.0-a1*(a1*mu+a1*mu*1.0/pow(a+s3,2.0)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+(la*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2)))*2.0+a*a11*k*2.0-a*a11*mu-a*a11*mu*1.0/pow(a+s1,2.0)+a*mu*1.0/pow(a+s1,3.0)*pow(a1+1.0,2.0)*2.0+la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*pow(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0),2.0)-a*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,3.0)*pow(a1+1.0,2.0)*2.0-(mu*((a1*a1)*(a*2.0+s1*2.0)+a11*(a+s1)*(a+s2)+a11*(a+s1)*(a+s3)+a11*(a+s2)*(a+s3)+a1*(a1+1.0)*(a*2.0+s2*2.0)+a1*(a1+1.0)*(a*2.0+s3*2.0)))/((a+s1)*(a+s2)*(a+s3))+a*a11*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)+(a1*mu*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2))+(a1*mu*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3))+(mu*1.0/pow(a+s1,2.0)*(a1+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3))-(a*la*1.0/pow(a+s1,2.0)*((a1*a1)*(a*2.0+s1*2.0)+a11*(a+s1)*(a+s2)+a11*(a+s1)*(a+s3)+a11*(a+s2)*(a+s3)+a1*(a1+1.0)*(a*2.0+s2*2.0)+a1*(a1+1.0)*(a*2.0+s3*2.0)))/((a+s2)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*((a1*a1)*(a*2.0+s1*2.0)+a11*(a+s1)*(a+s2)+a11*(a+s1)*(a+s3)+a11*(a+s2)*(a+s3)+a1*(a1+1.0)*(a*2.0+s2*2.0)+a1*(a1+1.0)*(a*2.0+s3*2.0)))/((a+s1)*(a+s2)*(a+s3))+(a*a1*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s2)+(a*a1*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s3)-(a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2))-(a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3))+(a*la*1.0/pow(a+s1,3.0)*(a1+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*3.0)/((a+s2)*(a+s3))-(la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a1+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3));
    dP_dF.x2222 = a2*(mu*(a2+1.0)+mu*1.0/pow(a+s2,2.0)*(a2+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a2+1.0)+(la*1.0/pow(a+s2,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3)))*-2.0-a22*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a22*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a22*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+(a2*a2)*k*2.0+mu*(pow(a2+1.0,2.0)*2.0+a22*(a*2.0+s1*2.0)+a22*(a*2.0+s2*2.0)+a22*(a*2.0+s3*2.0)+(a2*a2)*4.0)*(1.0/2.0)+a*(-a22*mu+(a2*a2)*mu*1.0/pow(a+s1,3.0)*2.0-a22*mu*1.0/pow(a+s1,2.0)+a22*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)-(a2*a2)*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,3.0)*2.0-(la*1.0/pow(a+s1,2.0)*((a2*a2)*(a*2.0+s2*2.0)+a22*(a+s1)*(a+s2)+a22*(a+s1)*(a+s3)+a22*(a+s2)*(a+s3)+a2*(a2+1.0)*(a*2.0+s1*2.0)+a2*(a2+1.0)*(a*2.0+s3*2.0)))/((a+s2)*(a+s3))+(a2*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s2)+(a2*la*1.0/pow(a+s1,3.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0))*3.0)/((a+s2)*(a+s3))+(la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a2+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s3))+a*(-a22*mu+(a2*a2)*mu*1.0/pow(a+s3,3.0)*2.0-a22*mu*1.0/pow(a+s3,2.0)+a22*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)-(a2*a2)*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,3.0)*2.0-(la*1.0/pow(a+s3,2.0)*((a2*a2)*(a*2.0+s2*2.0)+a22*(a+s1)*(a+s2)+a22*(a+s1)*(a+s3)+a22*(a+s2)*(a+s3)+a2*(a2+1.0)*(a*2.0+s1*2.0)+a2*(a2+1.0)*(a*2.0+s3*2.0)))/((a+s1)*(a+s2))+(a2*la*1.0/pow(a+s3,3.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0))*3.0)/((a+s1)*(a+s2))+(a2*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s2)+(la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a2+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s1))-a2*(a2*mu+a2*mu*1.0/pow(a+s1,2.0)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)+(la*1.0/pow(a+s1,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3)))*2.0-a2*(a2*mu+a2*mu*1.0/pow(a+s3,2.0)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+(la*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2)))*2.0+a*a22*k*2.0-a*a22*mu-a*a22*mu*1.0/pow(a+s2,2.0)+a*mu*1.0/pow(a+s2,3.0)*pow(a2+1.0,2.0)*2.0+la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*pow(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0),2.0)-a*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,3.0)*pow(a2+1.0,2.0)*2.0-(mu*((a2*a2)*(a*2.0+s2*2.0)+a22*(a+s1)*(a+s2)+a22*(a+s1)*(a+s3)+a22*(a+s2)*(a+s3)+a2*(a2+1.0)*(a*2.0+s1*2.0)+a2*(a2+1.0)*(a*2.0+s3*2.0)))/((a+s1)*(a+s2)*(a+s3))+a*a22*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+(a2*mu*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2))+(a2*mu*1.0/pow(a+s1,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3))+(mu*1.0/pow(a+s2,2.0)*(a2+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3))-(a*la*1.0/pow(a+s2,2.0)*((a2*a2)*(a*2.0+s2*2.0)+a22*(a+s1)*(a+s2)+a22*(a+s1)*(a+s3)+a22*(a+s2)*(a+s3)+a2*(a2+1.0)*(a*2.0+s1*2.0)+a2*(a2+1.0)*(a*2.0+s3*2.0)))/((a+s1)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*((a2*a2)*(a*2.0+s2*2.0)+a22*(a+s1)*(a+s2)+a22*(a+s1)*(a+s3)+a22*(a+s2)*(a+s3)+a2*(a2+1.0)*(a*2.0+s1*2.0)+a2*(a2+1.0)*(a*2.0+s3*2.0)))/((a+s1)*(a+s2)*(a+s3))+(a*a2*la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s1)+(a*a2*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s3)-(a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2))-(a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3))+(a*la*1.0/pow(a+s2,3.0)*(a2+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0))*3.0)/((a+s1)*(a+s3))-(la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a2+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3));
    dP_dF.x3333 = a3*(mu*(a3+1.0)+mu*1.0/pow(a+s3,2.0)*(a3+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a3+1.0)+(la*1.0/pow(a+s3,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2)))*-2.0-a33*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a33*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a33*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+(a3*a3)*k*2.0+mu*(pow(a3+1.0,2.0)*2.0+a33*(a*2.0+s1*2.0)+a33*(a*2.0+s2*2.0)+a33*(a*2.0+s3*2.0)+(a3*a3)*4.0)*(1.0/2.0)+a*(-a33*mu+(a3*a3)*mu*1.0/pow(a+s1,3.0)*2.0-a33*mu*1.0/pow(a+s1,2.0)+a33*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)-(a3*a3)*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,3.0)*2.0-(la*1.0/pow(a+s1,2.0)*((a3*a3)*(a*2.0+s3*2.0)+a33*(a+s1)*(a+s2)+a33*(a+s1)*(a+s3)+a33*(a+s2)*(a+s3)+a3*(a3+1.0)*(a*2.0+s1*2.0)+a3*(a3+1.0)*(a*2.0+s2*2.0)))/((a+s2)*(a+s3))+(a3*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/(a+s3)+(a3*la*1.0/pow(a+s1,3.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0))*3.0)/((a+s2)*(a+s3))+(la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a3+1.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/(a+s2))+a*(-a33*mu+(a3*a3)*mu*1.0/pow(a+s2,3.0)*2.0-a33*mu*1.0/pow(a+s2,2.0)+a33*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)-(a3*a3)*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,3.0)*2.0-(la*1.0/pow(a+s2,2.0)*((a3*a3)*(a*2.0+s3*2.0)+a33*(a+s1)*(a+s2)+a33*(a+s1)*(a+s3)+a33*(a+s2)*(a+s3)+a3*(a3+1.0)*(a*2.0+s1*2.0)+a3*(a3+1.0)*(a*2.0+s2*2.0)))/((a+s1)*(a+s3))+(a3*la*1.0/pow(a+s2,3.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0))*3.0)/((a+s1)*(a+s3))+(a3*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/(a+s3)+(la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a3+1.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/(a+s1))-a3*(a3*mu+a3*mu*1.0/pow(a+s1,2.0)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)+(la*1.0/pow(a+s1,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s2)*(a+s3)))*2.0-a3*(a3*mu+a3*mu*1.0/pow(a+s2,2.0)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+(la*1.0/pow(a+s2,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s3)))*2.0+a*a33*k*2.0-a*a33*mu-a*a33*mu*1.0/pow(a+s3,2.0)+a*mu*1.0/pow(a+s3,3.0)*pow(a3+1.0,2.0)*2.0+la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*pow(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0),2.0)-a*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,3.0)*pow(a3+1.0,2.0)*2.0-(mu*((a3*a3)*(a*2.0+s3*2.0)+a33*(a+s1)*(a+s2)+a33*(a+s1)*(a+s3)+a33*(a+s2)*(a+s3)+a3*(a3+1.0)*(a*2.0+s1*2.0)+a3*(a3+1.0)*(a*2.0+s2*2.0)))/((a+s1)*(a+s2)*(a+s3))+a*a33*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+(a3*mu*1.0/pow(a+s2,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s3))+(a3*mu*1.0/pow(a+s1,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s2)*(a+s3))+(mu*1.0/pow(a+s3,2.0)*(a3+1.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2))-(a*la*1.0/pow(a+s3,2.0)*((a3*a3)*(a*2.0+s3*2.0)+a33*(a+s1)*(a+s2)+a33*(a+s1)*(a+s3)+a33*(a+s2)*(a+s3)+a3*(a3+1.0)*(a*2.0+s1*2.0)+a3*(a3+1.0)*(a*2.0+s2*2.0)))/((a+s1)*(a+s2))+(la*log((a+s1)*(a+s2)*(a+s3))*((a3*a3)*(a*2.0+s3*2.0)+a33*(a+s1)*(a+s2)+a33*(a+s1)*(a+s3)+a33*(a+s2)*(a+s3)+a3*(a3+1.0)*(a*2.0+s1*2.0)+a3*(a3+1.0)*(a*2.0+s2*2.0)))/((a+s1)*(a+s2)*(a+s3))+(a*a3*la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/(a+s1)+(a*a3*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/(a+s2)-(a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s3))-(a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s2)*(a+s3))+(a*la*1.0/pow(a+s3,3.0)*(a3+1.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0))*3.0)/((a+s1)*(a+s2))-(la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a3+1.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2));
    dP_dF.x2211 = -a2*(mu*(a1+1.0)+mu*1.0/pow(a+s1,2.0)*(a1+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a1+1.0)+(la*1.0/pow(a+s1,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3)))-a1*(mu*(a2+1.0)+mu*1.0/pow(a+s2,2.0)*(a2+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a2+1.0)+(la*1.0/pow(a+s2,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3)))+mu*(a1*a2*2.0+a12*(a*2.0+s1*2.0)+a12*(a*2.0+s2*2.0)+a12*(a*2.0+s3*2.0)+a1*(a2*2.0+2.0)+a2*(a1*2.0+2.0))*(1.0/2.0)-a12*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a12*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a12*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+a*(-a12*mu-a12*mu*1.0/pow(a+s2,2.0)+a12*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+a1*mu*1.0/pow(a+s2,3.0)*(a2+1.0)*2.0-(la*1.0/pow(a+s2,2.0)*((a+s3)*(a1+1.0)*(a2+1.0)+a12*(a+s1)*(a+s2)+a12*(a+s1)*(a+s3)+a12*(a+s2)*(a+s3)+a1*a2*(a+s1)+a1*a2*(a+s2)+a1*a2*(a+s3)+a1*(a+s1)*(a2+1.0)+a2*(a+s2)*(a1+1.0)))/((a+s1)*(a+s3))-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,3.0)*(a2+1.0)*2.0+(a2*la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s1)+(a2*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s3)+(a1*la*1.0/pow(a+s2,3.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3))+(la*1.0/pow(a+s2,3.0)*(a2+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*2.0)/((a+s1)*(a+s3)))+a*(-a12*mu-a12*mu*1.0/pow(a+s3,2.0)+a1*a2*mu*1.0/pow(a+s3,3.0)*2.0+a12*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)-a1*a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,3.0)*2.0-(la*1.0/pow(a+s3,2.0)*((a+s3)*(a1+1.0)*(a2+1.0)+a12*(a+s1)*(a+s2)+a12*(a+s1)*(a+s3)+a12*(a+s2)*(a+s3)+a1*a2*(a+s1)+a1*a2*(a+s2)+a1*a2*(a+s3)+a1*(a+s1)*(a2+1.0)+a2*(a+s2)*(a1+1.0)))/((a+s1)*(a+s2))+(a2*la*1.0/pow(a+s3,3.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*2.0)/((a+s1)*(a+s2))+(a2*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s2)+(a1*la*1.0/pow(a+s3,3.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2))+(la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a2+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s1))-a2*(a1*mu+a1*mu*1.0/pow(a+s2,2.0)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+(la*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3)))-a2*(a1*mu+a1*mu*1.0/pow(a+s3,2.0)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+(la*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2)))-a1*(a2*mu+a2*mu*1.0/pow(a+s1,2.0)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)+(la*1.0/pow(a+s1,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3)))-a1*(a2*mu+a2*mu*1.0/pow(a+s3,2.0)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+(la*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2)))+a*a12*k*2.0+a1*a2*k*2.0-a*a12*mu-a*a12*mu*1.0/pow(a+s1,2.0)+a*a12*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)-(mu*((a+s3)*(a1+1.0)*(a2+1.0)+a12*(a+s1)*(a+s2)+a12*(a+s1)*(a+s3)+a12*(a+s2)*(a+s3)+a1*a2*(a+s1)+a1*a2*(a+s2)+a1*a2*(a+s3)+a1*(a+s1)*(a2+1.0)+a2*(a+s2)*(a1+1.0)))/((a+s1)*(a+s2)*(a+s3))+a*a2*mu*1.0/pow(a+s1,3.0)*(a1+1.0)*2.0+(a2*mu*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2))+(a2*mu*1.0/pow(a+s1,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3))-(a*la*1.0/pow(a+s1,2.0)*((a+s3)*(a1+1.0)*(a2+1.0)+a12*(a+s1)*(a+s2)+a12*(a+s1)*(a+s3)+a12*(a+s2)*(a+s3)+a1*a2*(a+s1)+a1*a2*(a+s2)+a1*a2*(a+s3)+a1*(a+s1)*(a2+1.0)+a2*(a+s2)*(a1+1.0)))/((a+s2)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*((a+s3)*(a1+1.0)*(a2+1.0)+a12*(a+s1)*(a+s2)+a12*(a+s1)*(a+s3)+a12*(a+s2)*(a+s3)+a1*a2*(a+s1)+a1*a2*(a+s2)+a1*a2*(a+s3)+a1*(a+s1)*(a2+1.0)+a2*(a+s2)*(a1+1.0)))/((a+s1)*(a+s2)*(a+s3))-a*a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,3.0)*(a1+1.0)*2.0+(mu*1.0/pow(a+s2,2.0)*(a2+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3))+la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0))+(a*a2*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s2)+(a*a2*la*1.0/pow(a+s1,3.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*2.0)/((a+s2)*(a+s3))-(a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2))-(a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3))+(a*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a2+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s3)+(a*la*1.0/pow(a+s1,3.0)*(a1+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3))-(la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a2+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3));
    dP_dF.x3311 = -a3*(mu*(a1+1.0)+mu*1.0/pow(a+s1,2.0)*(a1+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a1+1.0)+(la*1.0/pow(a+s1,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3)))-a1*(mu*(a3+1.0)+mu*1.0/pow(a+s3,2.0)*(a3+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a3+1.0)+(la*1.0/pow(a+s3,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2)))+mu*(a1*a3*2.0+a13*(a*2.0+s1*2.0)+a13*(a*2.0+s2*2.0)+a13*(a*2.0+s3*2.0)+a1*(a3*2.0+2.0)+a3*(a1*2.0+2.0))*(1.0/2.0)-a13*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a13*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a13*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+a*(-a13*mu-a13*mu*1.0/pow(a+s3,2.0)+a13*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+a1*mu*1.0/pow(a+s3,3.0)*(a3+1.0)*2.0-(la*1.0/pow(a+s3,2.0)*((a+s2)*(a1+1.0)*(a3+1.0)+a13*(a+s1)*(a+s2)+a13*(a+s1)*(a+s3)+a13*(a+s2)*(a+s3)+a1*a3*(a+s1)+a1*a3*(a+s2)+a1*a3*(a+s3)+a1*(a+s1)*(a3+1.0)+a3*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2))-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,3.0)*(a3+1.0)*2.0+(a3*la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s1)+(a3*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s2)+(a1*la*1.0/pow(a+s3,3.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2))+(la*1.0/pow(a+s3,3.0)*(a3+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*2.0)/((a+s1)*(a+s2)))+a*(-a13*mu-a13*mu*1.0/pow(a+s2,2.0)+a1*a3*mu*1.0/pow(a+s2,3.0)*2.0+a13*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)-a1*a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,3.0)*2.0-(la*1.0/pow(a+s2,2.0)*((a+s2)*(a1+1.0)*(a3+1.0)+a13*(a+s1)*(a+s2)+a13*(a+s1)*(a+s3)+a13*(a+s2)*(a+s3)+a1*a3*(a+s1)+a1*a3*(a+s2)+a1*a3*(a+s3)+a1*(a+s1)*(a3+1.0)+a3*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3))+(a3*la*1.0/pow(a+s2,3.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*2.0)/((a+s1)*(a+s3))+(a3*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s3)+(a1*la*1.0/pow(a+s2,3.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s3))+(la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a3+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s1))-a3*(a1*mu+a1*mu*1.0/pow(a+s2,2.0)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+(la*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3)))-a3*(a1*mu+a1*mu*1.0/pow(a+s3,2.0)-a1*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+(la*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2)))-a1*(a3*mu+a3*mu*1.0/pow(a+s1,2.0)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)+(la*1.0/pow(a+s1,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s2)*(a+s3)))-a1*(a3*mu+a3*mu*1.0/pow(a+s2,2.0)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+(la*1.0/pow(a+s2,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s3)))+a*a13*k*2.0+a1*a3*k*2.0-a*a13*mu-a*a13*mu*1.0/pow(a+s1,2.0)+a*a13*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)-(mu*((a+s2)*(a1+1.0)*(a3+1.0)+a13*(a+s1)*(a+s2)+a13*(a+s1)*(a+s3)+a13*(a+s2)*(a+s3)+a1*a3*(a+s1)+a1*a3*(a+s2)+a1*a3*(a+s3)+a1*(a+s1)*(a3+1.0)+a3*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2)*(a+s3))+a*a3*mu*1.0/pow(a+s1,3.0)*(a1+1.0)*2.0+(a3*mu*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3))+(a3*mu*1.0/pow(a+s1,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3))-(a*la*1.0/pow(a+s1,2.0)*((a+s2)*(a1+1.0)*(a3+1.0)+a13*(a+s1)*(a+s2)+a13*(a+s1)*(a+s3)+a13*(a+s2)*(a+s3)+a1*a3*(a+s1)+a1*a3*(a+s2)+a1*a3*(a+s3)+a1*(a+s1)*(a3+1.0)+a3*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*((a+s2)*(a1+1.0)*(a3+1.0)+a13*(a+s1)*(a+s2)+a13*(a+s1)*(a+s3)+a13*(a+s2)*(a+s3)+a1*a3*(a+s1)+a1*a3*(a+s2)+a1*a3*(a+s3)+a1*(a+s1)*(a3+1.0)+a3*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2)*(a+s3))-a*a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,3.0)*(a1+1.0)*2.0+(mu*1.0/pow(a+s3,2.0)*(a3+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2))+la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0))+(a*a3*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s3)+(a*a3*la*1.0/pow(a+s1,3.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0))*2.0)/((a+s2)*(a+s3))-(a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s3))-(a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s2)*(a+s3))+(a*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a3+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/(a+s2)+(a*la*1.0/pow(a+s1,3.0)*(a1+1.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s2)*(a+s3))-(la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a3+1.0)*(a1*(a+s1)*(a+s2)+a1*(a+s1)*(a+s3)+(a+s2)*(a+s3)*(a1+1.0)))/((a+s1)*(a+s2));
    dP_dF.x3322 = -a3*(mu*(a2+1.0)+mu*1.0/pow(a+s2,2.0)*(a2+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a2+1.0)+(la*1.0/pow(a+s2,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3)))-a2*(mu*(a3+1.0)+mu*1.0/pow(a+s3,2.0)*(a3+1.0)-la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a3+1.0)+(la*1.0/pow(a+s3,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2)))+mu*(a2*a3*2.0+a23*(a*2.0+s1*2.0)+a23*(a*2.0+s2*2.0)+a23*(a*2.0+s3*2.0)+a2*(a3*2.0+2.0)+a3*(a2*2.0+2.0))*(1.0/2.0)-a23*(mu*(a+s1)-mu/(a+s1)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s1))-a23*(mu*(a+s2)-mu/(a+s2)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s2))-a23*(mu*(a+s3)-mu/(a+s3)+(la*log((a+s1)*(a+s2)*(a+s3)))/(a+s3))+a*(-a23*mu-a23*mu*1.0/pow(a+s3,2.0)+a23*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+a2*mu*1.0/pow(a+s3,3.0)*(a3+1.0)*2.0-(la*1.0/pow(a+s3,2.0)*((a+s1)*(a2+1.0)*(a3+1.0)+a23*(a+s1)*(a+s2)+a23*(a+s1)*(a+s3)+a23*(a+s2)*(a+s3)+a2*a3*(a+s1)+a2*a3*(a+s2)+a2*a3*(a+s3)+a2*(a+s2)*(a3+1.0)+a3*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2))-a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,3.0)*(a3+1.0)*2.0+(a3*la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s1)+(a3*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s2)+(a2*la*1.0/pow(a+s3,3.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s2))+(la*1.0/pow(a+s3,3.0)*(a3+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0))*2.0)/((a+s1)*(a+s2)))+a*(-a23*mu-a23*mu*1.0/pow(a+s1,2.0)+a2*a3*mu*1.0/pow(a+s1,3.0)*2.0+a23*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)-a2*a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,3.0)*2.0-(la*1.0/pow(a+s1,2.0)*((a+s1)*(a2+1.0)*(a3+1.0)+a23*(a+s1)*(a+s2)+a23*(a+s1)*(a+s3)+a23*(a+s2)*(a+s3)+a2*a3*(a+s1)+a2*a3*(a+s2)+a2*a3*(a+s3)+a2*(a+s2)*(a3+1.0)+a3*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3))+(a3*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s3)+(a3*la*1.0/pow(a+s1,3.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0))*2.0)/((a+s2)*(a+s3))+(a2*la*1.0/pow(a+s1,3.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s2)*(a+s3))+(la*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0)*(a3+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s2))-a3*(a2*mu+a2*mu*1.0/pow(a+s1,2.0)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)+(la*1.0/pow(a+s1,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3)))-a3*(a2*mu+a2*mu*1.0/pow(a+s3,2.0)-a2*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)+(la*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2)))-a2*(a3*mu+a3*mu*1.0/pow(a+s1,2.0)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)+(la*1.0/pow(a+s1,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s2)*(a+s3)))-a2*(a3*mu+a3*mu*1.0/pow(a+s2,2.0)-a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)+(la*1.0/pow(a+s2,2.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s3)))+a*a23*k*2.0+a2*a3*k*2.0-a*a23*mu-a*a23*mu*1.0/pow(a+s2,2.0)+a*a23*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)-(mu*((a+s1)*(a2+1.0)*(a3+1.0)+a23*(a+s1)*(a+s2)+a23*(a+s1)*(a+s3)+a23*(a+s2)*(a+s3)+a2*a3*(a+s1)+a2*a3*(a+s2)+a2*a3*(a+s3)+a2*(a+s2)*(a3+1.0)+a3*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2)*(a+s3))+a*a3*mu*1.0/pow(a+s2,3.0)*(a2+1.0)*2.0+(a3*mu*1.0/pow(a+s2,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3))+(a3*mu*1.0/pow(a+s1,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3))-(a*la*1.0/pow(a+s2,2.0)*((a+s1)*(a2+1.0)*(a3+1.0)+a23*(a+s1)*(a+s2)+a23*(a+s1)*(a+s3)+a23*(a+s2)*(a+s3)+a2*a3*(a+s1)+a2*a3*(a+s2)+a2*a3*(a+s3)+a2*(a+s2)*(a3+1.0)+a3*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3))+(la*log((a+s1)*(a+s2)*(a+s3))*((a+s1)*(a2+1.0)*(a3+1.0)+a23*(a+s1)*(a+s2)+a23*(a+s1)*(a+s3)+a23*(a+s2)*(a+s3)+a2*a3*(a+s1)+a2*a3*(a+s2)+a2*a3*(a+s3)+a2*(a+s2)*(a3+1.0)+a3*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2)*(a+s3))-a*a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,3.0)*(a2+1.0)*2.0+(mu*1.0/pow(a+s3,2.0)*(a3+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2))+la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0))*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0))+(a*a3*la*1.0/pow(a+s2,3.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0))*2.0)/((a+s1)*(a+s3))+(a*a3*la*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s3)-(a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s2,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s3))-(a3*la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s1,2.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s2)*(a+s3))+(a*la*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a3+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/(a+s1)+(a*la*1.0/pow(a+s2,3.0)*(a2+1.0)*(a3*(a+s1)*(a+s3)+a3*(a+s2)*(a+s3)+(a+s1)*(a+s2)*(a3+1.0)))/((a+s1)*(a+s3))-(la*log((a+s1)*(a+s2)*(a+s3))*1.0/pow(a+s3,2.0)*(a3+1.0)*(a2*(a+s1)*(a+s2)+a2*(a+s2)*(a+s3)+(a+s1)*(a+s3)*(a2+1.0)))/((a+s1)*(a+s2));
    dP_dF.x2121 = ((a*a*a)*mu*s1*2.0-(a*a*a)*mu*s2*2.0+(a*a)*la*(s1*s1)-(a*a)*la*(s2*s2)+(a*a)*mu*(s1*s1)-(a*a)*mu*(s1*s1*s1*s1)-(a*a*a)*mu*(s1*s1*s1)*2.0-(a*a*a*a)*mu*(s1*s1)-(a*a)*mu*(s2*s2)+(a*a)*mu*(s2*s2*s2*s2)+(a*a*a)*mu*(s2*s2*s2)*2.0+(a*a*a*a)*mu*(s2*s2)+mu*(s1*s1)*(s2*s2*s2*s2)-mu*(s1*s1*s1*s1)*(s2*s2)+(a*a*a)*la*s1*2.0-(a*a*a)*la*s2*2.0-(a*a)*la*(s1*s1)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)+(a*a)*la*(s2*s2)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)-(a*a*a*a*a)*a1*k*s1*2.0+(a*a*a*a*a)*a2*k*s2*2.0+a*a1*la*(s1*s1*s1)+(a*a*a)*a1*la*s1*4.0-a*a2*la*(s2*s2*s2)-(a*a*a)*a2*la*s2*4.0+a*a1*mu*(s1*s1*s1)+(a*a*a)*a1*mu*s1*2.0+(a*a*a*a*a)*a1*mu*s1*3.0-a*a2*mu*(s2*s2*s2)-(a*a*a)*a2*mu*s2*2.0-(a*a*a*a*a)*a2*mu*s2*3.0+a*mu*s1*(s2*s2*s2*s2)*2.0-a*mu*(s1*s1*s1*s1)*s2*2.0-(a*a*a)*a1*k*(s1*s1*s1)*2.0-(a*a*a*a)*a1*k*(s1*s1)*4.0+(a*a*a)*a2*k*(s2*s2*s2)*2.0+(a*a*a*a)*a2*k*(s2*s2)*4.0+(a*a)*a1*la*(s1*s1)*4.0-(a*a)*a2*la*(s2*s2)*4.0+(a*a)*a1*mu*(s1*s1)*2.0+(a*a*a)*a1*mu*(s1*s1*s1)*3.0+(a*a*a*a)*a1*mu*(s1*s1)*6.0-(a*a)*a2*mu*(s2*s2)*2.0-(a*a*a)*a2*mu*(s2*s2*s2)*3.0-(a*a*a*a)*a2*mu*(s2*s2)*6.0+a*mu*(s1*s1)*(s2*s2*s2)*2.0-a*mu*(s1*s1*s1)*(s2*s2)*2.0+(a*a)*mu*s1*(s2*s2*s2)*4.0-(a*a)*mu*(s1*s1*s1)*s2*4.0+(a*a*a)*mu*s1*(s2*s2)*2.0-(a*a*a)*mu*(s1*s1)*s2*2.0-(a*a*a)*la*s1*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+(a*a*a)*la*s2*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0-a*a1*k*(s1*s1*s1)*(s2*s2)*2.0-(a*a)*a1*k*(s1*s1*s1)*s2*4.0-(a*a*a)*a1*k*s1*(s2*s2)*2.0-(a*a*a)*a1*k*(s1*s1)*s2*8.0+a*a2*k*(s1*s1)*(s2*s2*s2)*2.0+(a*a)*a2*k*s1*(s2*s2*s2)*4.0+(a*a*a)*a2*k*s1*(s2*s2)*8.0+(a*a*a)*a2*k*(s1*s1)*s2*2.0+a*a1*mu*(s1*s1*s1)*(s2*s2)*3.0+(a*a)*a1*mu*(s1*s1*s1)*s2*6.0+(a*a*a)*a1*mu*s1*(s2*s2)*3.0+(a*a*a)*a1*mu*(s1*s1)*s2*1.2E1-a*a2*mu*(s1*s1)*(s2*s2*s2)*3.0-(a*a)*a2*mu*s1*(s2*s2*s2)*6.0-(a*a*a)*a2*mu*s1*(s2*s2)*1.2E1-(a*a*a)*a2*mu*(s1*s1)*s2*3.0-a*a1*la*(s1*s1*s1)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)-(a*a*a)*a1*la*s1*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+a*a2*la*(s2*s2*s2)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)+(a*a*a)*a2*la*s2*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0-(a*a)*a1*k*(s1*s1)*(s2*s2)*4.0+(a*a)*a2*k*(s1*s1)*(s2*s2)*4.0+(a*a)*a1*mu*(s1*s1)*(s2*s2)*6.0-(a*a)*a2*mu*(s1*s1)*(s2*s2)*6.0-(a*a)*a1*la*(s1*s1)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+(a*a)*a2*la*(s2*s2)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0-(a*a*a*a)*a1*k*s1*s2*4.0+(a*a*a*a)*a2*k*s1*s2*4.0+a*a1*la*s1*(s2*s2)+a*a1*la*(s1*s1)*s2*2.0+(a*a)*a1*la*s1*s2*4.0-a*a2*la*s1*(s2*s2)*2.0-a*a2*la*(s1*s1)*s2-(a*a)*a2*la*s1*s2*4.0+a*a1*mu*s1*(s2*s2)+(a*a)*a1*mu*s1*s2*2.0+(a*a*a*a)*a1*mu*s1*s2*6.0-a*a2*mu*(s1*s1)*s2-(a*a)*a2*mu*s1*s2*2.0-(a*a*a*a)*a2*mu*s1*s2*6.0-a*a1*la*s1*(s2*s2)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)-(a*a)*a1*la*s1*s2*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+a*a2*la*(s1*s1)*s2*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)+(a*a)*a2*la*s1*s2*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0)/((s1*s1)*(s2*s2*s2*s2)-(s1*s1*s1*s1)*(s2*s2)-(a*a)*(s1*s1*s1*s1)-(a*a*a)*(s1*s1*s1)*2.0-(a*a*a*a)*(s1*s1)+(a*a)*(s2*s2*s2*s2)+(a*a*a)*(s2*s2*s2)*2.0+(a*a*a*a)*(s2*s2)+a*s1*(s2*s2*s2*s2)*2.0-a*(s1*s1*s1*s1)*s2*2.0+a*(s1*s1)*(s2*s2*s2)*2.0-a*(s1*s1*s1)*(s2*s2)*2.0+(a*a)*s1*(s2*s2*s2)*4.0-(a*a)*(s1*s1*s1)*s2*4.0+(a*a*a)*s1*(s2*s2)*2.0-(a*a*a)*(s1*s1)*s2*2.0)-((a*pow(a+s1,2.0)*pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*(a1*s1-a2*s2)*(la+mu-la*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3))+(a*la*(a+s1)*(a+s2)*(a*s1-a*s2+a1*(s1*s1)*2.0-a2*(s2*s2)*2.0+a*a1*s1*4.0-a*a2*s2*4.0+a1*s1*s2*2.0-a2*s1*s2*2.0))/(a+s3))*1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0))/((s1+s2)*(s1-s2));
    dP_dF.x3131 = ((a*a*a)*mu*s1*2.0-(a*a*a)*mu*s3*2.0+(a*a)*la*(s1*s1)-(a*a)*la*(s3*s3)+(a*a)*mu*(s1*s1)-(a*a)*mu*(s1*s1*s1*s1)-(a*a*a)*mu*(s1*s1*s1)*2.0-(a*a*a*a)*mu*(s1*s1)-(a*a)*mu*(s3*s3)+(a*a)*mu*(s3*s3*s3*s3)+(a*a*a)*mu*(s3*s3*s3)*2.0+(a*a*a*a)*mu*(s3*s3)+mu*(s1*s1)*(s3*s3*s3*s3)-mu*(s1*s1*s1*s1)*(s3*s3)+(a*a*a)*la*s1*2.0-(a*a*a)*la*s3*2.0-(a*a)*la*(s1*s1)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)+(a*a)*la*(s3*s3)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)-(a*a*a*a*a)*a1*k*s1*2.0+(a*a*a*a*a)*a3*k*s3*2.0+a*a1*la*(s1*s1*s1)+(a*a*a)*a1*la*s1*4.0-a*a3*la*(s3*s3*s3)-(a*a*a)*a3*la*s3*4.0+a*a1*mu*(s1*s1*s1)+(a*a*a)*a1*mu*s1*2.0+(a*a*a*a*a)*a1*mu*s1*3.0-a*a3*mu*(s3*s3*s3)-(a*a*a)*a3*mu*s3*2.0-(a*a*a*a*a)*a3*mu*s3*3.0+a*mu*s1*(s3*s3*s3*s3)*2.0-a*mu*(s1*s1*s1*s1)*s3*2.0-(a*a*a)*a1*k*(s1*s1*s1)*2.0-(a*a*a*a)*a1*k*(s1*s1)*4.0+(a*a*a)*a3*k*(s3*s3*s3)*2.0+(a*a*a*a)*a3*k*(s3*s3)*4.0+(a*a)*a1*la*(s1*s1)*4.0-(a*a)*a3*la*(s3*s3)*4.0+(a*a)*a1*mu*(s1*s1)*2.0+(a*a*a)*a1*mu*(s1*s1*s1)*3.0+(a*a*a*a)*a1*mu*(s1*s1)*6.0-(a*a)*a3*mu*(s3*s3)*2.0-(a*a*a)*a3*mu*(s3*s3*s3)*3.0-(a*a*a*a)*a3*mu*(s3*s3)*6.0+a*mu*(s1*s1)*(s3*s3*s3)*2.0-a*mu*(s1*s1*s1)*(s3*s3)*2.0+(a*a)*mu*s1*(s3*s3*s3)*4.0-(a*a)*mu*(s1*s1*s1)*s3*4.0+(a*a*a)*mu*s1*(s3*s3)*2.0-(a*a*a)*mu*(s1*s1)*s3*2.0-(a*a*a)*la*s1*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+(a*a*a)*la*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0-a*a1*k*(s1*s1*s1)*(s3*s3)*2.0-(a*a)*a1*k*(s1*s1*s1)*s3*4.0-(a*a*a)*a1*k*s1*(s3*s3)*2.0-(a*a*a)*a1*k*(s1*s1)*s3*8.0+a*a3*k*(s1*s1)*(s3*s3*s3)*2.0+(a*a)*a3*k*s1*(s3*s3*s3)*4.0+(a*a*a)*a3*k*s1*(s3*s3)*8.0+(a*a*a)*a3*k*(s1*s1)*s3*2.0+a*a1*mu*(s1*s1*s1)*(s3*s3)*3.0+(a*a)*a1*mu*(s1*s1*s1)*s3*6.0+(a*a*a)*a1*mu*s1*(s3*s3)*3.0+(a*a*a)*a1*mu*(s1*s1)*s3*1.2E1-a*a3*mu*(s1*s1)*(s3*s3*s3)*3.0-(a*a)*a3*mu*s1*(s3*s3*s3)*6.0-(a*a*a)*a3*mu*s1*(s3*s3)*1.2E1-(a*a*a)*a3*mu*(s1*s1)*s3*3.0-a*a1*la*(s1*s1*s1)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)-(a*a*a)*a1*la*s1*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+a*a3*la*(s3*s3*s3)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)+(a*a*a)*a3*la*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0-(a*a)*a1*k*(s1*s1)*(s3*s3)*4.0+(a*a)*a3*k*(s1*s1)*(s3*s3)*4.0+(a*a)*a1*mu*(s1*s1)*(s3*s3)*6.0-(a*a)*a3*mu*(s1*s1)*(s3*s3)*6.0-(a*a)*a1*la*(s1*s1)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+(a*a)*a3*la*(s3*s3)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0-(a*a*a*a)*a1*k*s1*s3*4.0+(a*a*a*a)*a3*k*s1*s3*4.0+a*a1*la*s1*(s3*s3)+a*a1*la*(s1*s1)*s3*2.0+(a*a)*a1*la*s1*s3*4.0-a*a3*la*s1*(s3*s3)*2.0-a*a3*la*(s1*s1)*s3-(a*a)*a3*la*s1*s3*4.0+a*a1*mu*s1*(s3*s3)+(a*a)*a1*mu*s1*s3*2.0+(a*a*a*a)*a1*mu*s1*s3*6.0-a*a3*mu*(s1*s1)*s3-(a*a)*a3*mu*s1*s3*2.0-(a*a*a*a)*a3*mu*s1*s3*6.0-a*a1*la*s1*(s3*s3)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)-(a*a)*a1*la*s1*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+a*a3*la*(s1*s1)*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)+(a*a)*a3*la*s1*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0)/((s1*s1)*(s3*s3*s3*s3)-(s1*s1*s1*s1)*(s3*s3)-(a*a)*(s1*s1*s1*s1)-(a*a*a)*(s1*s1*s1)*2.0-(a*a*a*a)*(s1*s1)+(a*a)*(s3*s3*s3*s3)+(a*a*a)*(s3*s3*s3)*2.0+(a*a*a*a)*(s3*s3)+a*s1*(s3*s3*s3*s3)*2.0-a*(s1*s1*s1*s1)*s3*2.0+a*(s1*s1)*(s3*s3*s3)*2.0-a*(s1*s1*s1)*(s3*s3)*2.0+(a*a)*s1*(s3*s3*s3)*4.0-(a*a)*(s1*s1*s1)*s3*4.0+(a*a*a)*s1*(s3*s3)*2.0-(a*a*a)*(s1*s1)*s3*2.0)-((a*pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*pow(a+s3,2.0)*(a1*s1-a3*s3)*(la+mu-la*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3))+(a*la*(a+s1)*(a+s3)*(a*s1-a*s3+a1*(s1*s1)*2.0-a3*(s3*s3)*2.0+a*a1*s1*4.0-a*a3*s3*4.0+a1*s1*s3*2.0-a3*s1*s3*2.0))/(a+s2))*1.0/pow(a+s1,2.0)*1.0/pow(a+s3,2.0))/((s1+s3)*(s1-s3));
    dP_dF.x3232 = ((a*a*a)*mu*s2*2.0-(a*a*a)*mu*s3*2.0+(a*a)*la*(s2*s2)-(a*a)*la*(s3*s3)+(a*a)*mu*(s2*s2)-(a*a)*mu*(s2*s2*s2*s2)-(a*a*a)*mu*(s2*s2*s2)*2.0-(a*a*a*a)*mu*(s2*s2)-(a*a)*mu*(s3*s3)+(a*a)*mu*(s3*s3*s3*s3)+(a*a*a)*mu*(s3*s3*s3)*2.0+(a*a*a*a)*mu*(s3*s3)+mu*(s2*s2)*(s3*s3*s3*s3)-mu*(s2*s2*s2*s2)*(s3*s3)+(a*a*a)*la*s2*2.0-(a*a*a)*la*s3*2.0-(a*a)*la*(s2*s2)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)+(a*a)*la*(s3*s3)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)-(a*a*a*a*a)*a2*k*s2*2.0+(a*a*a*a*a)*a3*k*s3*2.0+a*a2*la*(s2*s2*s2)+(a*a*a)*a2*la*s2*4.0-a*a3*la*(s3*s3*s3)-(a*a*a)*a3*la*s3*4.0+a*a2*mu*(s2*s2*s2)+(a*a*a)*a2*mu*s2*2.0+(a*a*a*a*a)*a2*mu*s2*3.0-a*a3*mu*(s3*s3*s3)-(a*a*a)*a3*mu*s3*2.0-(a*a*a*a*a)*a3*mu*s3*3.0+a*mu*s2*(s3*s3*s3*s3)*2.0-a*mu*(s2*s2*s2*s2)*s3*2.0-(a*a*a)*a2*k*(s2*s2*s2)*2.0-(a*a*a*a)*a2*k*(s2*s2)*4.0+(a*a*a)*a3*k*(s3*s3*s3)*2.0+(a*a*a*a)*a3*k*(s3*s3)*4.0+(a*a)*a2*la*(s2*s2)*4.0-(a*a)*a3*la*(s3*s3)*4.0+(a*a)*a2*mu*(s2*s2)*2.0+(a*a*a)*a2*mu*(s2*s2*s2)*3.0+(a*a*a*a)*a2*mu*(s2*s2)*6.0-(a*a)*a3*mu*(s3*s3)*2.0-(a*a*a)*a3*mu*(s3*s3*s3)*3.0-(a*a*a*a)*a3*mu*(s3*s3)*6.0+a*mu*(s2*s2)*(s3*s3*s3)*2.0-a*mu*(s2*s2*s2)*(s3*s3)*2.0+(a*a)*mu*s2*(s3*s3*s3)*4.0-(a*a)*mu*(s2*s2*s2)*s3*4.0+(a*a*a)*mu*s2*(s3*s3)*2.0-(a*a*a)*mu*(s2*s2)*s3*2.0-(a*a*a)*la*s2*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+(a*a*a)*la*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0-a*a2*k*(s2*s2*s2)*(s3*s3)*2.0-(a*a)*a2*k*(s2*s2*s2)*s3*4.0-(a*a*a)*a2*k*s2*(s3*s3)*2.0-(a*a*a)*a2*k*(s2*s2)*s3*8.0+a*a3*k*(s2*s2)*(s3*s3*s3)*2.0+(a*a)*a3*k*s2*(s3*s3*s3)*4.0+(a*a*a)*a3*k*s2*(s3*s3)*8.0+(a*a*a)*a3*k*(s2*s2)*s3*2.0+a*a2*mu*(s2*s2*s2)*(s3*s3)*3.0+(a*a)*a2*mu*(s2*s2*s2)*s3*6.0+(a*a*a)*a2*mu*s2*(s3*s3)*3.0+(a*a*a)*a2*mu*(s2*s2)*s3*1.2E1-a*a3*mu*(s2*s2)*(s3*s3*s3)*3.0-(a*a)*a3*mu*s2*(s3*s3*s3)*6.0-(a*a*a)*a3*mu*s2*(s3*s3)*1.2E1-(a*a*a)*a3*mu*(s2*s2)*s3*3.0-a*a2*la*(s2*s2*s2)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)-(a*a*a)*a2*la*s2*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+a*a3*la*(s3*s3*s3)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)+(a*a*a)*a3*la*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0-(a*a)*a2*k*(s2*s2)*(s3*s3)*4.0+(a*a)*a3*k*(s2*s2)*(s3*s3)*4.0+(a*a)*a2*mu*(s2*s2)*(s3*s3)*6.0-(a*a)*a3*mu*(s2*s2)*(s3*s3)*6.0-(a*a)*a2*la*(s2*s2)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+(a*a)*a3*la*(s3*s3)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0-(a*a*a*a)*a2*k*s2*s3*4.0+(a*a*a*a)*a3*k*s2*s3*4.0+a*a2*la*s2*(s3*s3)+a*a2*la*(s2*s2)*s3*2.0+(a*a)*a2*la*s2*s3*4.0-a*a3*la*s2*(s3*s3)*2.0-a*a3*la*(s2*s2)*s3-(a*a)*a3*la*s2*s3*4.0+a*a2*mu*s2*(s3*s3)+(a*a)*a2*mu*s2*s3*2.0+(a*a*a*a)*a2*mu*s2*s3*6.0-a*a3*mu*(s2*s2)*s3-(a*a)*a3*mu*s2*s3*2.0-(a*a*a*a)*a3*mu*s2*s3*6.0-a*a2*la*s2*(s3*s3)*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)-(a*a)*a2*la*s2*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0+a*a3*la*(s2*s2)*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)+(a*a)*a3*la*s2*s3*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3)*2.0)/((s2*s2)*(s3*s3*s3*s3)-(s2*s2*s2*s2)*(s3*s3)-(a*a)*(s2*s2*s2*s2)-(a*a*a)*(s2*s2*s2)*2.0-(a*a*a*a)*(s2*s2)+(a*a)*(s3*s3*s3*s3)+(a*a*a)*(s3*s3*s3)*2.0+(a*a*a*a)*(s3*s3)+a*s2*(s3*s3*s3*s3)*2.0-a*(s2*s2*s2*s2)*s3*2.0+a*(s2*s2)*(s3*s3*s3)*2.0-a*(s2*s2*s2)*(s3*s3)*2.0+(a*a)*s2*(s3*s3*s3)*4.0-(a*a)*(s2*s2*s2)*s3*4.0+(a*a*a)*s2*(s3*s3)*2.0-(a*a*a)*(s2*s2)*s3*2.0)-((a*1.0/pow(a+s1,2.0)*pow(a+s2,2.0)*pow(a+s3,2.0)*(a2*s2-a3*s3)*(la+mu-la*log((a*a)*s1+(a*a)*s2+(a*a)*s3+a*a*a+a*s1*s2+a*s1*s3+a*s2*s3+s1*s2*s3))+(a*la*(a+s2)*(a+s3)*(a*s2-a*s3+a2*(s2*s2)*2.0-a3*(s3*s3)*2.0+a*a2*s2*4.0-a*a3*s3*4.0+a2*s2*s3*2.0-a3*s2*s3*2.0))/(a+s1))*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0))/((s2+s3)*(s2-s3));
    dP_dF.x2121 = (1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*((a*a*a*a*a)*mu*s1*2.0-(a*a*a*a*a)*mu*s2*2.0+(a*a*a)*la*(s1*s1*s1)*2.0+(a*a*a*a)*la*(s1*s1)*5.0-(a*a*a)*la*(s2*s2*s2)*2.0-(a*a*a*a)*la*(s2*s2)*5.0+(a*a*a)*mu*(s1*s1*s1)*2.0+(a*a*a*a)*mu*(s1*s1)*4.0-(a*a*a)*mu*(s2*s2*s2)*2.0-(a*a*a*a)*mu*(s2*s2)*4.0+(a*a*a*a*a)*la*s1*3.0-(a*a*a*a*a)*la*s2*3.0+(a*a)*la*(s1*s1)*(s3*s3)*3.0-(a*a)*la*(s2*s2)*(s3*s3)*3.0+(a*a)*mu*(s1*s1)*(s3*s3)*4.0-(a*a)*mu*(s2*s2)*(s3*s3)*4.0+(a*a*a*a*a*a*a)*a1*k*s2*2.0-(a*a*a*a*a*a*a)*a2*k*s1*2.0-(a*a*a*a*a)*a1*la*s2*9.0+(a*a*a*a*a)*a2*la*s1*9.0-(a*a*a*a*a)*a1*mu*s2*3.0+(a*a*a*a*a)*a2*mu*s1*3.0-(a*a*a*a*a*a*a)*a1*mu*s2*3.0+(a*a*a*a*a*a*a)*a2*mu*s1*3.0-(a*a*a*a*a)*la*s1*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a*a)*la*s2*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a)*la*s1*s3*5.0-(a*a*a*a)*la*s2*s3*5.0+(a*a*a*a)*mu*s1*s3*4.0-(a*a*a*a)*mu*s2*s3*4.0+(a*a*a*a*a)*a1*k*(s2*s2*s2)*2.0-(a*a*a*a*a)*a2*k*(s1*s1*s1)*2.0+(a*a*a*a*a*a)*a1*k*(s2*s2)*4.0-(a*a*a*a*a*a)*a2*k*(s1*s1)*4.0-(a*a*a)*a1*la*(s2*s2*s2)*4.0+(a*a*a)*a2*la*(s1*s1*s1)*4.0-(a*a*a*a)*a1*la*(s2*s2)*1.2E1+(a*a*a*a)*a2*la*(s1*s1)*1.2E1-(a*a*a)*a1*mu*(s2*s2*s2)*2.0+(a*a*a)*a2*mu*(s1*s1*s1)*2.0-(a*a*a*a)*a1*mu*(s2*s2)*4.0+(a*a*a*a)*a2*mu*(s1*s1)*4.0-(a*a*a*a*a)*a1*mu*(s2*s2*s2)*3.0+(a*a*a*a*a)*a2*mu*(s1*s1*s1)*3.0-(a*a*a*a*a*a)*a1*mu*(s2*s2)*6.0+(a*a*a*a*a*a)*a2*mu*(s1*s1)*6.0-(a*a*a)*la*(s1*s1*s1)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a*a)*la*(s1*s1)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*la*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a)*la*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a)*la*s1*(s2*s2*s2)+(a*a)*la*(s1*s1*s1)*s2-(a*a*a)*la*s1*(s2*s2)*3.0+(a*a*a)*la*(s1*s1)*s2*3.0+a*la*(s1*s1*s1)*(s3*s3)+(a*a)*la*(s1*s1*s1)*s3*3.0+(a*a*a)*la*s1*(s3*s3)*2.0+(a*a*a)*la*(s1*s1)*s3*8.0-a*la*(s2*s2*s2)*(s3*s3)-(a*a)*la*(s2*s2*s2)*s3*3.0-(a*a*a)*la*s2*(s3*s3)*2.0-(a*a*a)*la*(s2*s2)*s3*8.0-(a*a)*mu*s1*(s2*s2*s2)+(a*a)*mu*(s1*s1*s1)*s2-(a*a*a)*mu*s1*(s2*s2)*2.0+(a*a*a)*mu*(s1*s1)*s2*2.0+a*mu*(s1*s1*s1)*(s3*s3)*2.0+(a*a)*mu*(s1*s1*s1)*s3*4.0+(a*a*a)*mu*s1*(s3*s3)*2.0+(a*a*a)*mu*(s1*s1)*s3*8.0-a*mu*(s2*s2*s2)*(s3*s3)*2.0-(a*a)*mu*(s2*s2*s2)*s3*4.0-(a*a*a)*mu*s2*(s3*s3)*2.0-(a*a*a)*mu*(s2*s2)*s3*8.0-mu*s1*(s2*s2*s2)*(s3*s3)+mu*(s1*s1*s1)*s2*(s3*s3)+(a*a*a)*a1*la*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a2*la*(s1*s1*s1)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a)*a1*la*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a*a)*a2*la*(s1*s1)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a)*a1*k*s1*(s2*s2*s2)*4.0+(a*a*a*a*a)*a1*k*s1*(s2*s2)*8.0+(a*a*a*a*a)*a1*k*(s1*s1)*s2*2.0-(a*a*a*a)*a2*k*(s1*s1*s1)*s2*4.0-(a*a*a*a*a)*a2*k*s1*(s2*s2)*2.0-(a*a*a*a*a)*a2*k*(s1*s1)*s2*8.0+(a*a*a*a)*a1*k*(s2*s2*s2)*s3*4.0-(a*a*a*a)*a2*k*(s1*s1*s1)*s3*4.0+(a*a*a*a*a)*a1*k*s2*(s3*s3)*2.0+(a*a*a*a*a)*a1*k*(s2*s2)*s3*8.0-(a*a*a*a*a)*a2*k*s1*(s3*s3)*2.0-(a*a*a*a*a)*a2*k*(s1*s1)*s3*8.0-a*a1*la*(s1*s1)*(s2*s2*s2)-(a*a)*a1*la*s1*(s2*s2*s2)*4.0-(a*a*a)*a1*la*s1*(s2*s2)*1.4E1-(a*a*a)*a1*la*(s1*s1)*s2*4.0+a*a2*la*(s1*s1*s1)*(s2*s2)+(a*a)*a2*la*(s1*s1*s1)*s2*4.0+(a*a*a)*a2*la*s1*(s2*s2)*4.0+(a*a*a)*a2*la*(s1*s1)*s2*1.4E1-a*a1*la*(s2*s2*s2)*(s3*s3)+a*a2*la*(s1*s1*s1)*(s3*s3)-(a*a)*a1*la*(s2*s2*s2)*s3*4.0+(a*a)*a2*la*(s1*s1*s1)*s3*4.0-(a*a*a)*a1*la*s2*(s3*s3)*4.0-(a*a*a)*a1*la*(s2*s2)*s3*1.4E1+(a*a*a)*a2*la*s1*(s3*s3)*4.0+(a*a*a)*a2*la*(s1*s1)*s3*1.4E1-a*a1*mu*(s1*s1)*(s2*s2*s2)-(a*a)*a1*mu*s1*(s2*s2*s2)*2.0-(a*a*a)*a1*mu*s1*(s2*s2)*4.0-(a*a*a)*a1*mu*(s1*s1)*s2*2.0-(a*a*a*a)*a1*mu*s1*(s2*s2*s2)*6.0-(a*a*a*a*a)*a1*mu*s1*(s2*s2)*1.2E1-(a*a*a*a*a)*a1*mu*(s1*s1)*s2*3.0+a*a2*mu*(s1*s1*s1)*(s2*s2)+(a*a)*a2*mu*(s1*s1*s1)*s2*2.0+(a*a*a)*a2*mu*s1*(s2*s2)*2.0+(a*a*a)*a2*mu*(s1*s1)*s2*4.0+(a*a*a*a)*a2*mu*(s1*s1*s1)*s2*6.0+(a*a*a*a*a)*a2*mu*s1*(s2*s2)*3.0+(a*a*a*a*a)*a2*mu*(s1*s1)*s2*1.2E1-a*a1*mu*(s2*s2*s2)*(s3*s3)+a*a2*mu*(s1*s1*s1)*(s3*s3)-(a*a)*a1*mu*(s2*s2*s2)*s3*2.0+(a*a)*a2*mu*(s1*s1*s1)*s3*2.0-(a*a*a)*a1*mu*s2*(s3*s3)*2.0-(a*a*a)*a1*mu*(s2*s2)*s3*4.0+(a*a*a)*a2*mu*s1*(s3*s3)*2.0+(a*a*a)*a2*mu*(s1*s1)*s3*4.0-(a*a*a*a)*a1*mu*(s2*s2*s2)*s3*6.0+(a*a*a*a)*a2*mu*(s1*s1*s1)*s3*6.0-(a*a*a*a*a)*a1*mu*s2*(s3*s3)*3.0-(a*a*a*a*a)*a1*mu*(s2*s2)*s3*1.2E1+(a*a*a*a*a)*a2*mu*s1*(s3*s3)*3.0+(a*a*a*a*a)*a2*mu*(s1*s1)*s3*1.2E1+(a*a)*la*s1*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))-(a*a)*la*(s1*s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))+(a*a*a)*la*s1*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*la*(s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*2.0-a*la*(s1*s1*s1)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*la*(s1*s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a)*la*s1*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*la*(s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*8.0+a*la*(s2*s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*la*(s2*s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*la*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*la*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*8.0-a*la*s1*(s2*s2)*(s3*s3)+a*la*(s1*s1)*s2*(s3*s3)-(a*a)*la*s1*(s2*s2)*s3*4.0+(a*a)*la*(s1*s1)*s2*s3*4.0-a*mu*s1*(s2*s2)*(s3*s3)*2.0+a*mu*(s1*s1)*s2*(s3*s3)*2.0-(a*a)*mu*s1*(s2*s2)*s3*4.0+(a*a)*mu*(s1*s1)*s2*s3*4.0+la*s1*(s2*s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))-la*(s1*s1*s1)*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))+(a*a*a)*a1*k*(s1*s1)*(s2*s2*s2)*2.0+(a*a*a*a)*a1*k*(s1*s1)*(s2*s2)*4.0-(a*a*a)*a2*k*(s1*s1*s1)*(s2*s2)*2.0-(a*a*a*a)*a2*k*(s1*s1)*(s2*s2)*4.0+(a*a*a)*a1*k*(s2*s2*s2)*(s3*s3)*2.0-(a*a*a)*a2*k*(s1*s1*s1)*(s3*s3)*2.0+(a*a*a*a)*a1*k*(s2*s2)*(s3*s3)*4.0-(a*a*a*a)*a2*k*(s1*s1)*(s3*s3)*4.0-(a*a)*a1*la*(s1*s1)*(s2*s2)*4.0+(a*a)*a2*la*(s1*s1)*(s2*s2)*4.0-(a*a)*a1*la*(s2*s2)*(s3*s3)*4.0+(a*a)*a2*la*(s1*s1)*(s3*s3)*4.0-(a*a)*a1*mu*(s1*s1)*(s2*s2)*2.0-(a*a*a)*a1*mu*(s1*s1)*(s2*s2*s2)*3.0-(a*a*a*a)*a1*mu*(s1*s1)*(s2*s2)*6.0+(a*a)*a2*mu*(s1*s1)*(s2*s2)*2.0+(a*a*a)*a2*mu*(s1*s1*s1)*(s2*s2)*3.0+(a*a*a*a)*a2*mu*(s1*s1)*(s2*s2)*6.0-(a*a)*a1*mu*(s2*s2)*(s3*s3)*2.0+(a*a)*a2*mu*(s1*s1)*(s3*s3)*2.0-(a*a*a)*a1*mu*(s2*s2*s2)*(s3*s3)*3.0+(a*a*a)*a2*mu*(s1*s1*s1)*(s3*s3)*3.0-(a*a*a*a)*a1*mu*(s2*s2)*(s3*s3)*6.0+(a*a*a*a)*a2*mu*(s1*s1)*(s3*s3)*6.0-(a*a)*la*(s1*s1)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a)*la*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a*a)*a1*la*s2*log((a+s1)*(a+s2)*(a+s3))*3.0-(a*a*a*a*a)*a2*la*s1*log((a+s1)*(a+s2)*(a+s3))*3.0+(a*a*a*a*a*a)*a1*k*s1*s2*4.0-(a*a*a*a*a*a)*a2*k*s1*s2*4.0+(a*a*a*a*a*a)*a1*k*s2*s3*4.0-(a*a*a*a*a*a)*a2*k*s1*s3*4.0-(a*a*a*a)*a1*la*s1*s2*1.2E1+(a*a*a*a)*a2*la*s1*s2*1.2E1-(a*a*a*a)*a1*la*s2*s3*1.2E1+(a*a*a*a)*a2*la*s1*s3*1.2E1-(a*a*a*a)*a1*mu*s1*s2*4.0-(a*a*a*a*a*a)*a1*mu*s1*s2*6.0+(a*a*a*a)*a2*mu*s1*s2*4.0+(a*a*a*a*a*a)*a2*mu*s1*s2*6.0-(a*a*a*a)*a1*mu*s2*s3*4.0+(a*a*a*a)*a2*mu*s1*s3*4.0-(a*a*a*a*a*a)*a1*mu*s2*s3*6.0+(a*a*a*a*a*a)*a2*mu*s1*s3*6.0-(a*a*a*a)*la*s1*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a)*la*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-a*la*s1*(s2*s2*s2)*s3+a*la*(s1*s1*s1)*s2*s3-a*mu*s1*(s2*s2*s2)*s3*2.0+a*mu*(s1*s1*s1)*s2*s3*2.0+a*a1*k*(s1*s1)*(s2*s2*s2)*(s3*s3)*2.0+(a*a)*a1*k*s1*(s2*s2*s2)*(s3*s3)*4.0+(a*a)*a1*k*(s1*s1)*(s2*s2*s2)*s3*4.0+(a*a*a)*a1*k*s1*(s2*s2)*(s3*s3)*8.0+(a*a*a)*a1*k*(s1*s1)*s2*(s3*s3)*2.0+(a*a*a)*a1*k*(s1*s1)*(s2*s2)*s3*8.0-a*a2*k*(s1*s1*s1)*(s2*s2)*(s3*s3)*2.0-(a*a)*a2*k*(s1*s1*s1)*s2*(s3*s3)*4.0-(a*a)*a2*k*(s1*s1*s1)*(s2*s2)*s3*4.0-(a*a*a)*a2*k*s1*(s2*s2)*(s3*s3)*2.0-(a*a*a)*a2*k*(s1*s1)*s2*(s3*s3)*8.0-(a*a*a)*a2*k*(s1*s1)*(s2*s2)*s3*8.0-a*a1*mu*(s1*s1)*(s2*s2*s2)*(s3*s3)*3.0-(a*a)*a1*mu*s1*(s2*s2*s2)*(s3*s3)*6.0-(a*a)*a1*mu*(s1*s1)*(s2*s2*s2)*s3*6.0-(a*a*a)*a1*mu*s1*(s2*s2)*(s3*s3)*1.2E1-(a*a*a)*a1*mu*(s1*s1)*s2*(s3*s3)*3.0-(a*a*a)*a1*mu*(s1*s1)*(s2*s2)*s3*1.2E1+a*a2*mu*(s1*s1*s1)*(s2*s2)*(s3*s3)*3.0+(a*a)*a2*mu*(s1*s1*s1)*s2*(s3*s3)*6.0+(a*a)*a2*mu*(s1*s1*s1)*(s2*s2)*s3*6.0+(a*a*a)*a2*mu*s1*(s2*s2)*(s3*s3)*3.0+(a*a*a)*a2*mu*(s1*s1)*s2*(s3*s3)*1.2E1+(a*a*a)*a2*mu*(s1*s1)*(s2*s2)*s3*1.2E1+(a*a*a*a)*a1*la*s1*s2*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a*a)*a2*la*s1*s2*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a)*a1*la*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a*a)*a2*la*s1*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a*a)*a1*k*s1*s2*s3*8.0-(a*a*a*a*a)*a2*k*s1*s2*s3*8.0-a*a1*la*s1*(s2*s2*s2)*s3*2.0-(a*a*a)*a1*la*s1*s2*s3*1.4E1+a*a2*la*(s1*s1*s1)*s2*s3*2.0+(a*a*a)*a2*la*s1*s2*s3*1.4E1-(a*a*a)*a1*mu*s1*s2*s3*4.0-(a*a*a*a*a)*a1*mu*s1*s2*s3*1.2E1+(a*a*a)*a2*mu*s1*s2*s3*4.0+(a*a*a*a*a)*a2*mu*s1*s2*s3*1.2E1+a*la*s1*(s2*s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-a*la*(s1*s1*s1)*s2*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*a1*k*(s1*s1)*(s2*s2)*(s3*s3)*4.0-(a*a)*a2*k*(s1*s1)*(s2*s2)*(s3*s3)*4.0-(a*a)*a1*mu*(s1*s1)*(s2*s2)*(s3*s3)*6.0+(a*a)*a2*mu*(s1*s1)*(s2*s2)*(s3*s3)*6.0+a*a1*la*(s1*s1)*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))+(a*a)*a1*la*s1*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a1*la*s1*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*a1*la*(s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*2.0-a*a2*la*(s1*s1*s1)*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))-(a*a)*a2*la*(s1*s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a2*la*s1*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a2*la*(s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*4.0+a*a1*la*(s2*s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))-a*a2*la*(s1*s1*s1)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))+(a*a)*a1*la*(s2*s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a2*la*(s1*s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a1*la*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a1*la*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a)*a2*la*s1*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a2*la*(s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*a1*k*s1*(s2*s2*s2)*s3*8.0+(a*a*a*a)*a1*k*s1*s2*(s3*s3)*4.0+(a*a*a*a)*a1*k*s1*(s2*s2)*s3*1.6E1+(a*a*a*a)*a1*k*(s1*s1)*s2*s3*4.0-(a*a*a)*a2*k*(s1*s1*s1)*s2*s3*8.0-(a*a*a*a)*a2*k*s1*s2*(s3*s3)*4.0-(a*a*a*a)*a2*k*s1*(s2*s2)*s3*4.0-(a*a*a*a)*a2*k*(s1*s1)*s2*s3*1.6E1-a*a1*la*s1*(s2*s2)*(s3*s3)*2.0-a*a1*la*(s1*s1)*s2*(s3*s3)-a*a1*la*(s1*s1)*(s2*s2)*s3*2.0-(a*a)*a1*la*s1*s2*(s3*s3)*4.0-(a*a)*a1*la*s1*(s2*s2)*s3*1.2E1-(a*a)*a1*la*(s1*s1)*s2*s3*4.0+a*a2*la*s1*(s2*s2)*(s3*s3)+a*a2*la*(s1*s1)*s2*(s3*s3)*2.0+a*a2*la*(s1*s1)*(s2*s2)*s3*2.0+(a*a)*a2*la*s1*s2*(s3*s3)*4.0+(a*a)*a2*la*s1*(s2*s2)*s3*4.0+(a*a)*a2*la*(s1*s1)*s2*s3*1.2E1-a*a1*mu*(s1*s1)*s2*(s3*s3)-(a*a)*a1*mu*s1*s2*(s3*s3)*2.0-(a*a)*a1*mu*(s1*s1)*s2*s3*2.0-(a*a*a)*a1*mu*s1*(s2*s2*s2)*s3*1.2E1-(a*a*a*a)*a1*mu*s1*s2*(s3*s3)*6.0-(a*a*a*a)*a1*mu*s1*(s2*s2)*s3*2.4E1-(a*a*a*a)*a1*mu*(s1*s1)*s2*s3*6.0+a*a2*mu*s1*(s2*s2)*(s3*s3)+(a*a)*a2*mu*s1*s2*(s3*s3)*2.0+(a*a)*a2*mu*s1*(s2*s2)*s3*2.0+(a*a*a)*a2*mu*(s1*s1*s1)*s2*s3*1.2E1+(a*a*a*a)*a2*mu*s1*s2*(s3*s3)*6.0+(a*a*a*a)*a2*mu*s1*(s2*s2)*s3*6.0+(a*a*a*a)*a2*mu*(s1*s1)*s2*s3*2.4E1+a*la*s1*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-a*la*(s1*s1)*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*la*s1*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a)*la*(s1*s1)*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a)*a1*la*(s1*s1)*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a2*la*(s1*s1)*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*a1*la*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a2*la*(s1*s1)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+a*a1*la*(s1*s1)*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))+(a*a)*a1*la*s1*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*a1*la*(s1*s1)*s2*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-a*a2*la*s1*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))-(a*a)*a2*la*s1*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a2*la*s1*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a1*la*s1*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a)*a2*la*s1*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0))/(s1*s1-s2*s2);
    dP_dF.x3113 = (1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*((a*a*a*a*a)*mu*s1*2.0-(a*a*a*a*a)*mu*s3*2.0+(a*a*a)*la*(s1*s1*s1)*2.0+(a*a*a*a)*la*(s1*s1)*5.0-(a*a*a)*la*(s3*s3*s3)*2.0-(a*a*a*a)*la*(s3*s3)*5.0+(a*a*a)*mu*(s1*s1*s1)*2.0+(a*a*a*a)*mu*(s1*s1)*4.0-(a*a*a)*mu*(s3*s3*s3)*2.0-(a*a*a*a)*mu*(s3*s3)*4.0+(a*a*a*a*a)*la*s1*3.0-(a*a*a*a*a)*la*s3*3.0+(a*a)*la*(s1*s1)*(s2*s2)*3.0-(a*a)*la*(s2*s2)*(s3*s3)*3.0+(a*a)*mu*(s1*s1)*(s2*s2)*4.0-(a*a)*mu*(s2*s2)*(s3*s3)*4.0+(a*a*a*a*a*a*a)*a1*k*s3*2.0-(a*a*a*a*a*a*a)*a3*k*s1*2.0-(a*a*a*a*a)*a1*la*s3*9.0+(a*a*a*a*a)*a3*la*s1*9.0-(a*a*a*a*a)*a1*mu*s3*3.0+(a*a*a*a*a)*a3*mu*s1*3.0-(a*a*a*a*a*a*a)*a1*mu*s3*3.0+(a*a*a*a*a*a*a)*a3*mu*s1*3.0-(a*a*a*a*a)*la*s1*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a*a)*la*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a)*la*s1*s2*5.0-(a*a*a*a)*la*s2*s3*5.0+(a*a*a*a)*mu*s1*s2*4.0-(a*a*a*a)*mu*s2*s3*4.0+(a*a*a*a*a)*a1*k*(s3*s3*s3)*2.0-(a*a*a*a*a)*a3*k*(s1*s1*s1)*2.0+(a*a*a*a*a*a)*a1*k*(s3*s3)*4.0-(a*a*a*a*a*a)*a3*k*(s1*s1)*4.0-(a*a*a)*a1*la*(s3*s3*s3)*4.0+(a*a*a)*a3*la*(s1*s1*s1)*4.0-(a*a*a*a)*a1*la*(s3*s3)*1.2E1+(a*a*a*a)*a3*la*(s1*s1)*1.2E1-(a*a*a)*a1*mu*(s3*s3*s3)*2.0+(a*a*a)*a3*mu*(s1*s1*s1)*2.0-(a*a*a*a)*a1*mu*(s3*s3)*4.0+(a*a*a*a)*a3*mu*(s1*s1)*4.0-(a*a*a*a*a)*a1*mu*(s3*s3*s3)*3.0+(a*a*a*a*a)*a3*mu*(s1*s1*s1)*3.0-(a*a*a*a*a*a)*a1*mu*(s3*s3)*6.0+(a*a*a*a*a*a)*a3*mu*(s1*s1)*6.0-(a*a*a)*la*(s1*s1*s1)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a*a)*la*(s1*s1)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*la*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a)*la*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+a*la*(s1*s1*s1)*(s2*s2)+(a*a)*la*(s1*s1*s1)*s2*3.0+(a*a*a)*la*s1*(s2*s2)*2.0+(a*a*a)*la*(s1*s1)*s2*8.0-(a*a)*la*s1*(s3*s3*s3)+(a*a)*la*(s1*s1*s1)*s3-(a*a*a)*la*s1*(s3*s3)*3.0+(a*a*a)*la*(s1*s1)*s3*3.0-a*la*(s2*s2)*(s3*s3*s3)-(a*a)*la*s2*(s3*s3*s3)*3.0-(a*a*a)*la*s2*(s3*s3)*8.0-(a*a*a)*la*(s2*s2)*s3*2.0+a*mu*(s1*s1*s1)*(s2*s2)*2.0+(a*a)*mu*(s1*s1*s1)*s2*4.0+(a*a*a)*mu*s1*(s2*s2)*2.0+(a*a*a)*mu*(s1*s1)*s2*8.0-(a*a)*mu*s1*(s3*s3*s3)+(a*a)*mu*(s1*s1*s1)*s3-(a*a*a)*mu*s1*(s3*s3)*2.0+(a*a*a)*mu*(s1*s1)*s3*2.0-a*mu*(s2*s2)*(s3*s3*s3)*2.0-(a*a)*mu*s2*(s3*s3*s3)*4.0-(a*a*a)*mu*s2*(s3*s3)*8.0-(a*a*a)*mu*(s2*s2)*s3*2.0-mu*s1*(s2*s2)*(s3*s3*s3)+mu*(s1*s1*s1)*(s2*s2)*s3+(a*a*a)*a1*la*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a3*la*(s1*s1*s1)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a)*a1*la*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a*a)*a3*la*(s1*s1)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a)*a1*k*s1*(s3*s3*s3)*4.0+(a*a*a*a*a)*a1*k*s1*(s3*s3)*8.0+(a*a*a*a*a)*a1*k*(s1*s1)*s3*2.0+(a*a*a*a)*a1*k*s2*(s3*s3*s3)*4.0-(a*a*a*a)*a3*k*(s1*s1*s1)*s2*4.0+(a*a*a*a*a)*a1*k*s2*(s3*s3)*8.0+(a*a*a*a*a)*a1*k*(s2*s2)*s3*2.0-(a*a*a*a*a)*a3*k*s1*(s2*s2)*2.0-(a*a*a*a*a)*a3*k*(s1*s1)*s2*8.0-(a*a*a*a)*a3*k*(s1*s1*s1)*s3*4.0-(a*a*a*a*a)*a3*k*s1*(s3*s3)*2.0-(a*a*a*a*a)*a3*k*(s1*s1)*s3*8.0-a*a1*la*(s1*s1)*(s3*s3*s3)-(a*a)*a1*la*s1*(s3*s3*s3)*4.0-(a*a*a)*a1*la*s1*(s3*s3)*1.4E1-(a*a*a)*a1*la*(s1*s1)*s3*4.0-a*a1*la*(s2*s2)*(s3*s3*s3)+a*a3*la*(s1*s1*s1)*(s2*s2)-(a*a)*a1*la*s2*(s3*s3*s3)*4.0+(a*a)*a3*la*(s1*s1*s1)*s2*4.0-(a*a*a)*a1*la*s2*(s3*s3)*1.4E1-(a*a*a)*a1*la*(s2*s2)*s3*4.0+(a*a*a)*a3*la*s1*(s2*s2)*4.0+(a*a*a)*a3*la*(s1*s1)*s2*1.4E1+a*a3*la*(s1*s1*s1)*(s3*s3)+(a*a)*a3*la*(s1*s1*s1)*s3*4.0+(a*a*a)*a3*la*s1*(s3*s3)*4.0+(a*a*a)*a3*la*(s1*s1)*s3*1.4E1-a*a1*mu*(s1*s1)*(s3*s3*s3)-(a*a)*a1*mu*s1*(s3*s3*s3)*2.0-(a*a*a)*a1*mu*s1*(s3*s3)*4.0-(a*a*a)*a1*mu*(s1*s1)*s3*2.0-(a*a*a*a)*a1*mu*s1*(s3*s3*s3)*6.0-(a*a*a*a*a)*a1*mu*s1*(s3*s3)*1.2E1-(a*a*a*a*a)*a1*mu*(s1*s1)*s3*3.0-a*a1*mu*(s2*s2)*(s3*s3*s3)+a*a3*mu*(s1*s1*s1)*(s2*s2)-(a*a)*a1*mu*s2*(s3*s3*s3)*2.0+(a*a)*a3*mu*(s1*s1*s1)*s2*2.0-(a*a*a)*a1*mu*s2*(s3*s3)*4.0-(a*a*a)*a1*mu*(s2*s2)*s3*2.0+(a*a*a)*a3*mu*s1*(s2*s2)*2.0+(a*a*a)*a3*mu*(s1*s1)*s2*4.0-(a*a*a*a)*a1*mu*s2*(s3*s3*s3)*6.0+(a*a*a*a)*a3*mu*(s1*s1*s1)*s2*6.0-(a*a*a*a*a)*a1*mu*s2*(s3*s3)*1.2E1-(a*a*a*a*a)*a1*mu*(s2*s2)*s3*3.0+(a*a*a*a*a)*a3*mu*s1*(s2*s2)*3.0+(a*a*a*a*a)*a3*mu*(s1*s1)*s2*1.2E1+a*a3*mu*(s1*s1*s1)*(s3*s3)+(a*a)*a3*mu*(s1*s1*s1)*s3*2.0+(a*a*a)*a3*mu*s1*(s3*s3)*2.0+(a*a*a)*a3*mu*(s1*s1)*s3*4.0+(a*a*a*a)*a3*mu*(s1*s1*s1)*s3*6.0+(a*a*a*a*a)*a3*mu*s1*(s3*s3)*3.0+(a*a*a*a*a)*a3*mu*(s1*s1)*s3*1.2E1-a*la*(s1*s1*s1)*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*la*(s1*s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a)*la*s1*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*la*(s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*8.0+(a*a)*la*s1*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))-(a*a)*la*(s1*s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))+(a*a*a)*la*s1*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*la*(s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+a*la*(s2*s2)*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*la*s2*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*la*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*8.0+(a*a*a)*la*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-a*la*s1*(s2*s2)*(s3*s3)+a*la*(s1*s1)*(s2*s2)*s3-(a*a)*la*s1*s2*(s3*s3)*4.0+(a*a)*la*(s1*s1)*s2*s3*4.0-a*mu*s1*(s2*s2)*(s3*s3)*2.0+a*mu*(s1*s1)*(s2*s2)*s3*2.0-(a*a)*mu*s1*s2*(s3*s3)*4.0+(a*a)*mu*(s1*s1)*s2*s3*4.0+la*s1*(s2*s2)*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))-la*(s1*s1*s1)*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))+(a*a*a)*a1*k*(s1*s1)*(s3*s3*s3)*2.0+(a*a*a*a)*a1*k*(s1*s1)*(s3*s3)*4.0+(a*a*a)*a1*k*(s2*s2)*(s3*s3*s3)*2.0-(a*a*a)*a3*k*(s1*s1*s1)*(s2*s2)*2.0+(a*a*a*a)*a1*k*(s2*s2)*(s3*s3)*4.0-(a*a*a*a)*a3*k*(s1*s1)*(s2*s2)*4.0-(a*a*a)*a3*k*(s1*s1*s1)*(s3*s3)*2.0-(a*a*a*a)*a3*k*(s1*s1)*(s3*s3)*4.0-(a*a)*a1*la*(s1*s1)*(s3*s3)*4.0-(a*a)*a1*la*(s2*s2)*(s3*s3)*4.0+(a*a)*a3*la*(s1*s1)*(s2*s2)*4.0+(a*a)*a3*la*(s1*s1)*(s3*s3)*4.0-(a*a)*a1*mu*(s1*s1)*(s3*s3)*2.0-(a*a*a)*a1*mu*(s1*s1)*(s3*s3*s3)*3.0-(a*a*a*a)*a1*mu*(s1*s1)*(s3*s3)*6.0-(a*a)*a1*mu*(s2*s2)*(s3*s3)*2.0+(a*a)*a3*mu*(s1*s1)*(s2*s2)*2.0-(a*a*a)*a1*mu*(s2*s2)*(s3*s3*s3)*3.0+(a*a*a)*a3*mu*(s1*s1*s1)*(s2*s2)*3.0-(a*a*a*a)*a1*mu*(s2*s2)*(s3*s3)*6.0+(a*a*a*a)*a3*mu*(s1*s1)*(s2*s2)*6.0+(a*a)*a3*mu*(s1*s1)*(s3*s3)*2.0+(a*a*a)*a3*mu*(s1*s1*s1)*(s3*s3)*3.0+(a*a*a*a)*a3*mu*(s1*s1)*(s3*s3)*6.0-(a*a)*la*(s1*s1)*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a)*la*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a*a)*a1*la*s3*log((a+s1)*(a+s2)*(a+s3))*3.0-(a*a*a*a*a)*a3*la*s1*log((a+s1)*(a+s2)*(a+s3))*3.0+(a*a*a*a*a*a)*a1*k*s1*s3*4.0+(a*a*a*a*a*a)*a1*k*s2*s3*4.0-(a*a*a*a*a*a)*a3*k*s1*s2*4.0-(a*a*a*a*a*a)*a3*k*s1*s3*4.0-(a*a*a*a)*a1*la*s1*s3*1.2E1-(a*a*a*a)*a1*la*s2*s3*1.2E1+(a*a*a*a)*a3*la*s1*s2*1.2E1+(a*a*a*a)*a3*la*s1*s3*1.2E1-(a*a*a*a)*a1*mu*s1*s3*4.0-(a*a*a*a*a*a)*a1*mu*s1*s3*6.0-(a*a*a*a)*a1*mu*s2*s3*4.0+(a*a*a*a)*a3*mu*s1*s2*4.0-(a*a*a*a*a*a)*a1*mu*s2*s3*6.0+(a*a*a*a*a*a)*a3*mu*s1*s2*6.0+(a*a*a*a)*a3*mu*s1*s3*4.0+(a*a*a*a*a*a)*a3*mu*s1*s3*6.0-(a*a*a*a)*la*s1*s2*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a)*la*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-a*la*s1*s2*(s3*s3*s3)+a*la*(s1*s1*s1)*s2*s3-a*mu*s1*s2*(s3*s3*s3)*2.0+a*mu*(s1*s1*s1)*s2*s3*2.0+a*a1*k*(s1*s1)*(s2*s2)*(s3*s3*s3)*2.0+(a*a)*a1*k*s1*(s2*s2)*(s3*s3*s3)*4.0+(a*a)*a1*k*(s1*s1)*s2*(s3*s3*s3)*4.0+(a*a*a)*a1*k*s1*(s2*s2)*(s3*s3)*8.0+(a*a*a)*a1*k*(s1*s1)*s2*(s3*s3)*8.0+(a*a*a)*a1*k*(s1*s1)*(s2*s2)*s3*2.0-a*a3*k*(s1*s1*s1)*(s2*s2)*(s3*s3)*2.0-(a*a)*a3*k*(s1*s1*s1)*s2*(s3*s3)*4.0-(a*a)*a3*k*(s1*s1*s1)*(s2*s2)*s3*4.0-(a*a*a)*a3*k*s1*(s2*s2)*(s3*s3)*2.0-(a*a*a)*a3*k*(s1*s1)*s2*(s3*s3)*8.0-(a*a*a)*a3*k*(s1*s1)*(s2*s2)*s3*8.0-a*a1*mu*(s1*s1)*(s2*s2)*(s3*s3*s3)*3.0-(a*a)*a1*mu*s1*(s2*s2)*(s3*s3*s3)*6.0-(a*a)*a1*mu*(s1*s1)*s2*(s3*s3*s3)*6.0-(a*a*a)*a1*mu*s1*(s2*s2)*(s3*s3)*1.2E1-(a*a*a)*a1*mu*(s1*s1)*s2*(s3*s3)*1.2E1-(a*a*a)*a1*mu*(s1*s1)*(s2*s2)*s3*3.0+a*a3*mu*(s1*s1*s1)*(s2*s2)*(s3*s3)*3.0+(a*a)*a3*mu*(s1*s1*s1)*s2*(s3*s3)*6.0+(a*a)*a3*mu*(s1*s1*s1)*(s2*s2)*s3*6.0+(a*a*a)*a3*mu*s1*(s2*s2)*(s3*s3)*3.0+(a*a*a)*a3*mu*(s1*s1)*s2*(s3*s3)*1.2E1+(a*a*a)*a3*mu*(s1*s1)*(s2*s2)*s3*1.2E1+(a*a*a*a)*a1*la*s1*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a)*a1*la*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a*a)*a3*la*s1*s2*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a*a)*a3*la*s1*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a*a)*a1*k*s1*s2*s3*8.0-(a*a*a*a*a)*a3*k*s1*s2*s3*8.0-a*a1*la*s1*s2*(s3*s3*s3)*2.0-(a*a*a)*a1*la*s1*s2*s3*1.4E1+a*a3*la*(s1*s1*s1)*s2*s3*2.0+(a*a*a)*a3*la*s1*s2*s3*1.4E1-(a*a*a)*a1*mu*s1*s2*s3*4.0-(a*a*a*a*a)*a1*mu*s1*s2*s3*1.2E1+(a*a*a)*a3*mu*s1*s2*s3*4.0+(a*a*a*a*a)*a3*mu*s1*s2*s3*1.2E1+a*la*s1*s2*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-a*la*(s1*s1*s1)*s2*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*a1*k*(s1*s1)*(s2*s2)*(s3*s3)*4.0-(a*a)*a3*k*(s1*s1)*(s2*s2)*(s3*s3)*4.0-(a*a)*a1*mu*(s1*s1)*(s2*s2)*(s3*s3)*6.0+(a*a)*a3*mu*(s1*s1)*(s2*s2)*(s3*s3)*6.0+a*a1*la*(s1*s1)*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))+(a*a)*a1*la*s1*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a1*la*s1*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*a1*la*(s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+a*a1*la*(s2*s2)*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))-a*a3*la*(s1*s1*s1)*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))+(a*a)*a1*la*s2*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a3*la*(s1*s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a1*la*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*a1*la*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a3*la*s1*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a3*la*(s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*4.0-a*a3*la*(s1*s1*s1)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))-(a*a)*a3*la*(s1*s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a3*la*s1*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a3*la*(s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*a1*k*s1*s2*(s3*s3*s3)*8.0+(a*a*a*a)*a1*k*s1*s2*(s3*s3)*1.6E1+(a*a*a*a)*a1*k*s1*(s2*s2)*s3*4.0+(a*a*a*a)*a1*k*(s1*s1)*s2*s3*4.0-(a*a*a)*a3*k*(s1*s1*s1)*s2*s3*8.0-(a*a*a*a)*a3*k*s1*s2*(s3*s3)*4.0-(a*a*a*a)*a3*k*s1*(s2*s2)*s3*4.0-(a*a*a*a)*a3*k*(s1*s1)*s2*s3*1.6E1-a*a1*la*s1*(s2*s2)*(s3*s3)*2.0-a*a1*la*(s1*s1)*s2*(s3*s3)*2.0-a*a1*la*(s1*s1)*(s2*s2)*s3-(a*a)*a1*la*s1*s2*(s3*s3)*1.2E1-(a*a)*a1*la*s1*(s2*s2)*s3*4.0-(a*a)*a1*la*(s1*s1)*s2*s3*4.0+a*a3*la*s1*(s2*s2)*(s3*s3)+a*a3*la*(s1*s1)*s2*(s3*s3)*2.0+a*a3*la*(s1*s1)*(s2*s2)*s3*2.0+(a*a)*a3*la*s1*s2*(s3*s3)*4.0+(a*a)*a3*la*s1*(s2*s2)*s3*4.0+(a*a)*a3*la*(s1*s1)*s2*s3*1.2E1-a*a1*mu*(s1*s1)*(s2*s2)*s3-(a*a)*a1*mu*s1*(s2*s2)*s3*2.0-(a*a)*a1*mu*(s1*s1)*s2*s3*2.0-(a*a*a)*a1*mu*s1*s2*(s3*s3*s3)*1.2E1-(a*a*a*a)*a1*mu*s1*s2*(s3*s3)*2.4E1-(a*a*a*a)*a1*mu*s1*(s2*s2)*s3*6.0-(a*a*a*a)*a1*mu*(s1*s1)*s2*s3*6.0+a*a3*mu*s1*(s2*s2)*(s3*s3)+(a*a)*a3*mu*s1*s2*(s3*s3)*2.0+(a*a)*a3*mu*s1*(s2*s2)*s3*2.0+(a*a*a)*a3*mu*(s1*s1*s1)*s2*s3*1.2E1+(a*a*a*a)*a3*mu*s1*s2*(s3*s3)*6.0+(a*a*a*a)*a3*mu*s1*(s2*s2)*s3*6.0+(a*a*a*a)*a3*mu*(s1*s1)*s2*s3*2.4E1+a*la*s1*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-a*la*(s1*s1)*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*la*s1*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a)*la*(s1*s1)*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a)*a1*la*(s1*s1)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*a1*la*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a3*la*(s1*s1)*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a3*la*(s1*s1)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+a*a1*la*(s1*s1)*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))+(a*a)*a1*la*s1*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*a1*la*(s1*s1)*s2*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-a*a3*la*s1*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))-(a*a)*a3*la*s1*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a3*la*s1*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a1*la*s1*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a)*a3*la*s1*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0))/(s1*s1-s3*s3);
    dP_dF.x3223 = (1.0/pow(a+s1,2.0)*1.0/pow(a+s2,2.0)*1.0/pow(a+s3,2.0)*((a*a*a*a*a)*mu*s2*2.0-(a*a*a*a*a)*mu*s3*2.0+(a*a*a)*la*(s2*s2*s2)*2.0+(a*a*a*a)*la*(s2*s2)*5.0-(a*a*a)*la*(s3*s3*s3)*2.0-(a*a*a*a)*la*(s3*s3)*5.0+(a*a*a)*mu*(s2*s2*s2)*2.0+(a*a*a*a)*mu*(s2*s2)*4.0-(a*a*a)*mu*(s3*s3*s3)*2.0-(a*a*a*a)*mu*(s3*s3)*4.0+(a*a*a*a*a)*la*s2*3.0-(a*a*a*a*a)*la*s3*3.0+(a*a)*la*(s1*s1)*(s2*s2)*3.0-(a*a)*la*(s1*s1)*(s3*s3)*3.0+(a*a)*mu*(s1*s1)*(s2*s2)*4.0-(a*a)*mu*(s1*s1)*(s3*s3)*4.0+(a*a*a*a*a*a*a)*a2*k*s3*2.0-(a*a*a*a*a*a*a)*a3*k*s2*2.0-(a*a*a*a*a)*a2*la*s3*9.0+(a*a*a*a*a)*a3*la*s2*9.0-(a*a*a*a*a)*a2*mu*s3*3.0+(a*a*a*a*a)*a3*mu*s2*3.0-(a*a*a*a*a*a*a)*a2*mu*s3*3.0+(a*a*a*a*a*a*a)*a3*mu*s2*3.0-(a*a*a*a*a)*la*s2*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a*a)*la*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a)*la*s1*s2*5.0-(a*a*a*a)*la*s1*s3*5.0+(a*a*a*a)*mu*s1*s2*4.0-(a*a*a*a)*mu*s1*s3*4.0+(a*a*a*a*a)*a2*k*(s3*s3*s3)*2.0-(a*a*a*a*a)*a3*k*(s2*s2*s2)*2.0+(a*a*a*a*a*a)*a2*k*(s3*s3)*4.0-(a*a*a*a*a*a)*a3*k*(s2*s2)*4.0-(a*a*a)*a2*la*(s3*s3*s3)*4.0+(a*a*a)*a3*la*(s2*s2*s2)*4.0-(a*a*a*a)*a2*la*(s3*s3)*1.2E1+(a*a*a*a)*a3*la*(s2*s2)*1.2E1-(a*a*a)*a2*mu*(s3*s3*s3)*2.0+(a*a*a)*a3*mu*(s2*s2*s2)*2.0-(a*a*a*a)*a2*mu*(s3*s3)*4.0+(a*a*a*a)*a3*mu*(s2*s2)*4.0-(a*a*a*a*a)*a2*mu*(s3*s3*s3)*3.0+(a*a*a*a*a)*a3*mu*(s2*s2*s2)*3.0-(a*a*a*a*a*a)*a2*mu*(s3*s3)*6.0+(a*a*a*a*a*a)*a3*mu*(s2*s2)*6.0-(a*a*a)*la*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a*a)*la*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*la*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a)*la*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+a*la*(s1*s1)*(s2*s2*s2)+(a*a)*la*s1*(s2*s2*s2)*3.0+(a*a*a)*la*s1*(s2*s2)*8.0+(a*a*a)*la*(s1*s1)*s2*2.0-a*la*(s1*s1)*(s3*s3*s3)-(a*a)*la*s1*(s3*s3*s3)*3.0-(a*a*a)*la*s1*(s3*s3)*8.0-(a*a*a)*la*(s1*s1)*s3*2.0-(a*a)*la*s2*(s3*s3*s3)+(a*a)*la*(s2*s2*s2)*s3-(a*a*a)*la*s2*(s3*s3)*3.0+(a*a*a)*la*(s2*s2)*s3*3.0+a*mu*(s1*s1)*(s2*s2*s2)*2.0+(a*a)*mu*s1*(s2*s2*s2)*4.0+(a*a*a)*mu*s1*(s2*s2)*8.0+(a*a*a)*mu*(s1*s1)*s2*2.0-a*mu*(s1*s1)*(s3*s3*s3)*2.0-(a*a)*mu*s1*(s3*s3*s3)*4.0-(a*a*a)*mu*s1*(s3*s3)*8.0-(a*a*a)*mu*(s1*s1)*s3*2.0-(a*a)*mu*s2*(s3*s3*s3)+(a*a)*mu*(s2*s2*s2)*s3-(a*a*a)*mu*s2*(s3*s3)*2.0+(a*a*a)*mu*(s2*s2)*s3*2.0-mu*(s1*s1)*s2*(s3*s3*s3)+mu*(s1*s1)*(s2*s2*s2)*s3+(a*a*a)*a2*la*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a3*la*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a*a)*a2*la*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a*a)*a3*la*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a)*a2*k*s1*(s3*s3*s3)*4.0-(a*a*a*a)*a3*k*s1*(s2*s2*s2)*4.0+(a*a*a*a*a)*a2*k*s1*(s3*s3)*8.0+(a*a*a*a*a)*a2*k*(s1*s1)*s3*2.0-(a*a*a*a*a)*a3*k*s1*(s2*s2)*8.0-(a*a*a*a*a)*a3*k*(s1*s1)*s2*2.0+(a*a*a*a)*a2*k*s2*(s3*s3*s3)*4.0+(a*a*a*a*a)*a2*k*s2*(s3*s3)*8.0+(a*a*a*a*a)*a2*k*(s2*s2)*s3*2.0-(a*a*a*a)*a3*k*(s2*s2*s2)*s3*4.0-(a*a*a*a*a)*a3*k*s2*(s3*s3)*2.0-(a*a*a*a*a)*a3*k*(s2*s2)*s3*8.0-a*a2*la*(s1*s1)*(s3*s3*s3)+a*a3*la*(s1*s1)*(s2*s2*s2)-(a*a)*a2*la*s1*(s3*s3*s3)*4.0+(a*a)*a3*la*s1*(s2*s2*s2)*4.0-(a*a*a)*a2*la*s1*(s3*s3)*1.4E1-(a*a*a)*a2*la*(s1*s1)*s3*4.0+(a*a*a)*a3*la*s1*(s2*s2)*1.4E1+(a*a*a)*a3*la*(s1*s1)*s2*4.0-a*a2*la*(s2*s2)*(s3*s3*s3)-(a*a)*a2*la*s2*(s3*s3*s3)*4.0-(a*a*a)*a2*la*s2*(s3*s3)*1.4E1-(a*a*a)*a2*la*(s2*s2)*s3*4.0+a*a3*la*(s2*s2*s2)*(s3*s3)+(a*a)*a3*la*(s2*s2*s2)*s3*4.0+(a*a*a)*a3*la*s2*(s3*s3)*4.0+(a*a*a)*a3*la*(s2*s2)*s3*1.4E1-a*a2*mu*(s1*s1)*(s3*s3*s3)+a*a3*mu*(s1*s1)*(s2*s2*s2)-(a*a)*a2*mu*s1*(s3*s3*s3)*2.0+(a*a)*a3*mu*s1*(s2*s2*s2)*2.0-(a*a*a)*a2*mu*s1*(s3*s3)*4.0-(a*a*a)*a2*mu*(s1*s1)*s3*2.0+(a*a*a)*a3*mu*s1*(s2*s2)*4.0+(a*a*a)*a3*mu*(s1*s1)*s2*2.0-(a*a*a*a)*a2*mu*s1*(s3*s3*s3)*6.0+(a*a*a*a)*a3*mu*s1*(s2*s2*s2)*6.0-(a*a*a*a*a)*a2*mu*s1*(s3*s3)*1.2E1-(a*a*a*a*a)*a2*mu*(s1*s1)*s3*3.0+(a*a*a*a*a)*a3*mu*s1*(s2*s2)*1.2E1+(a*a*a*a*a)*a3*mu*(s1*s1)*s2*3.0-a*a2*mu*(s2*s2)*(s3*s3*s3)-(a*a)*a2*mu*s2*(s3*s3*s3)*2.0-(a*a*a)*a2*mu*s2*(s3*s3)*4.0-(a*a*a)*a2*mu*(s2*s2)*s3*2.0-(a*a*a*a)*a2*mu*s2*(s3*s3*s3)*6.0-(a*a*a*a*a)*a2*mu*s2*(s3*s3)*1.2E1-(a*a*a*a*a)*a2*mu*(s2*s2)*s3*3.0+a*a3*mu*(s2*s2*s2)*(s3*s3)+(a*a)*a3*mu*(s2*s2*s2)*s3*2.0+(a*a*a)*a3*mu*s2*(s3*s3)*2.0+(a*a*a)*a3*mu*(s2*s2)*s3*4.0+(a*a*a*a)*a3*mu*(s2*s2*s2)*s3*6.0+(a*a*a*a*a)*a3*mu*s2*(s3*s3)*3.0+(a*a*a*a*a)*a3*mu*(s2*s2)*s3*1.2E1-a*la*(s1*s1)*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*la*s1*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a)*la*s1*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*8.0-(a*a*a)*la*(s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*2.0+a*la*(s1*s1)*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*la*s1*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*la*s1*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*8.0+(a*a*a)*la*(s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*la*s2*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))-(a*a)*la*(s2*s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))+(a*a*a)*la*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*la*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-a*la*(s1*s1)*s2*(s3*s3)+a*la*(s1*s1)*(s2*s2)*s3-(a*a)*la*s1*s2*(s3*s3)*4.0+(a*a)*la*s1*(s2*s2)*s3*4.0-a*mu*(s1*s1)*s2*(s3*s3)*2.0+a*mu*(s1*s1)*(s2*s2)*s3*2.0-(a*a)*mu*s1*s2*(s3*s3)*4.0+(a*a)*mu*s1*(s2*s2)*s3*4.0+la*(s1*s1)*s2*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))-la*(s1*s1)*(s2*s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))+(a*a*a)*a2*k*(s1*s1)*(s3*s3*s3)*2.0-(a*a*a)*a3*k*(s1*s1)*(s2*s2*s2)*2.0+(a*a*a*a)*a2*k*(s1*s1)*(s3*s3)*4.0-(a*a*a*a)*a3*k*(s1*s1)*(s2*s2)*4.0+(a*a*a)*a2*k*(s2*s2)*(s3*s3*s3)*2.0+(a*a*a*a)*a2*k*(s2*s2)*(s3*s3)*4.0-(a*a*a)*a3*k*(s2*s2*s2)*(s3*s3)*2.0-(a*a*a*a)*a3*k*(s2*s2)*(s3*s3)*4.0-(a*a)*a2*la*(s1*s1)*(s3*s3)*4.0+(a*a)*a3*la*(s1*s1)*(s2*s2)*4.0-(a*a)*a2*la*(s2*s2)*(s3*s3)*4.0+(a*a)*a3*la*(s2*s2)*(s3*s3)*4.0-(a*a)*a2*mu*(s1*s1)*(s3*s3)*2.0+(a*a)*a3*mu*(s1*s1)*(s2*s2)*2.0-(a*a*a)*a2*mu*(s1*s1)*(s3*s3*s3)*3.0+(a*a*a)*a3*mu*(s1*s1)*(s2*s2*s2)*3.0-(a*a*a*a)*a2*mu*(s1*s1)*(s3*s3)*6.0+(a*a*a*a)*a3*mu*(s1*s1)*(s2*s2)*6.0-(a*a)*a2*mu*(s2*s2)*(s3*s3)*2.0-(a*a*a)*a2*mu*(s2*s2)*(s3*s3*s3)*3.0-(a*a*a*a)*a2*mu*(s2*s2)*(s3*s3)*6.0+(a*a)*a3*mu*(s2*s2)*(s3*s3)*2.0+(a*a*a)*a3*mu*(s2*s2*s2)*(s3*s3)*3.0+(a*a*a*a)*a3*mu*(s2*s2)*(s3*s3)*6.0-(a*a)*la*(s1*s1)*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a)*la*(s1*s1)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a*a)*a2*la*s3*log((a+s1)*(a+s2)*(a+s3))*3.0-(a*a*a*a*a)*a3*la*s2*log((a+s1)*(a+s2)*(a+s3))*3.0+(a*a*a*a*a*a)*a2*k*s1*s3*4.0-(a*a*a*a*a*a)*a3*k*s1*s2*4.0+(a*a*a*a*a*a)*a2*k*s2*s3*4.0-(a*a*a*a*a*a)*a3*k*s2*s3*4.0-(a*a*a*a)*a2*la*s1*s3*1.2E1+(a*a*a*a)*a3*la*s1*s2*1.2E1-(a*a*a*a)*a2*la*s2*s3*1.2E1+(a*a*a*a)*a3*la*s2*s3*1.2E1-(a*a*a*a)*a2*mu*s1*s3*4.0+(a*a*a*a)*a3*mu*s1*s2*4.0-(a*a*a*a*a*a)*a2*mu*s1*s3*6.0+(a*a*a*a*a*a)*a3*mu*s1*s2*6.0-(a*a*a*a)*a2*mu*s2*s3*4.0-(a*a*a*a*a*a)*a2*mu*s2*s3*6.0+(a*a*a*a)*a3*mu*s2*s3*4.0+(a*a*a*a*a*a)*a3*mu*s2*s3*6.0-(a*a*a*a)*la*s1*s2*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a)*la*s1*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-a*la*s1*s2*(s3*s3*s3)+a*la*s1*(s2*s2*s2)*s3-a*mu*s1*s2*(s3*s3*s3)*2.0+a*mu*s1*(s2*s2*s2)*s3*2.0+a*a2*k*(s1*s1)*(s2*s2)*(s3*s3*s3)*2.0+(a*a)*a2*k*s1*(s2*s2)*(s3*s3*s3)*4.0+(a*a)*a2*k*(s1*s1)*s2*(s3*s3*s3)*4.0+(a*a*a)*a2*k*s1*(s2*s2)*(s3*s3)*8.0+(a*a*a)*a2*k*(s1*s1)*s2*(s3*s3)*8.0+(a*a*a)*a2*k*(s1*s1)*(s2*s2)*s3*2.0-a*a3*k*(s1*s1)*(s2*s2*s2)*(s3*s3)*2.0-(a*a)*a3*k*s1*(s2*s2*s2)*(s3*s3)*4.0-(a*a)*a3*k*(s1*s1)*(s2*s2*s2)*s3*4.0-(a*a*a)*a3*k*s1*(s2*s2)*(s3*s3)*8.0-(a*a*a)*a3*k*(s1*s1)*s2*(s3*s3)*2.0-(a*a*a)*a3*k*(s1*s1)*(s2*s2)*s3*8.0-a*a2*mu*(s1*s1)*(s2*s2)*(s3*s3*s3)*3.0-(a*a)*a2*mu*s1*(s2*s2)*(s3*s3*s3)*6.0-(a*a)*a2*mu*(s1*s1)*s2*(s3*s3*s3)*6.0-(a*a*a)*a2*mu*s1*(s2*s2)*(s3*s3)*1.2E1-(a*a*a)*a2*mu*(s1*s1)*s2*(s3*s3)*1.2E1-(a*a*a)*a2*mu*(s1*s1)*(s2*s2)*s3*3.0+a*a3*mu*(s1*s1)*(s2*s2*s2)*(s3*s3)*3.0+(a*a)*a3*mu*s1*(s2*s2*s2)*(s3*s3)*6.0+(a*a)*a3*mu*(s1*s1)*(s2*s2*s2)*s3*6.0+(a*a*a)*a3*mu*s1*(s2*s2)*(s3*s3)*1.2E1+(a*a*a)*a3*mu*(s1*s1)*s2*(s3*s3)*3.0+(a*a*a)*a3*mu*(s1*s1)*(s2*s2)*s3*1.2E1+(a*a*a*a)*a2*la*s1*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a*a)*a3*la*s1*s2*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a)*a2*la*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a*a)*a3*la*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a*a*a)*a2*k*s1*s2*s3*8.0-(a*a*a*a*a)*a3*k*s1*s2*s3*8.0-a*a2*la*s1*s2*(s3*s3*s3)*2.0-(a*a*a)*a2*la*s1*s2*s3*1.4E1+a*a3*la*s1*(s2*s2*s2)*s3*2.0+(a*a*a)*a3*la*s1*s2*s3*1.4E1-(a*a*a)*a2*mu*s1*s2*s3*4.0-(a*a*a*a*a)*a2*mu*s1*s2*s3*1.2E1+(a*a*a)*a3*mu*s1*s2*s3*4.0+(a*a*a*a*a)*a3*mu*s1*s2*s3*1.2E1+a*la*s1*s2*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-a*la*s1*(s2*s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*a2*k*(s1*s1)*(s2*s2)*(s3*s3)*4.0-(a*a)*a3*k*(s1*s1)*(s2*s2)*(s3*s3)*4.0-(a*a)*a2*mu*(s1*s1)*(s2*s2)*(s3*s3)*6.0+(a*a)*a3*mu*(s1*s1)*(s2*s2)*(s3*s3)*6.0+a*a2*la*(s1*s1)*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))-a*a3*la*(s1*s1)*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))+(a*a)*a2*la*s1*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a3*la*s1*(s2*s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a2*la*s1*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*a2*la*(s1*s1)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a3*la*s1*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a)*a3*la*(s1*s1)*s2*log((a+s1)*(a+s2)*(a+s3))*2.0+a*a2*la*(s2*s2)*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))+(a*a)*a2*la*s2*(s3*s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a2*la*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*a2*la*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-a*a3*la*(s2*s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))-(a*a)*a3*la*(s2*s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a3*la*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a*a)*a3*la*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a*a)*a2*k*s1*s2*(s3*s3*s3)*8.0+(a*a*a*a)*a2*k*s1*s2*(s3*s3)*1.6E1+(a*a*a*a)*a2*k*s1*(s2*s2)*s3*4.0+(a*a*a*a)*a2*k*(s1*s1)*s2*s3*4.0-(a*a*a)*a3*k*s1*(s2*s2*s2)*s3*8.0-(a*a*a*a)*a3*k*s1*s2*(s3*s3)*4.0-(a*a*a*a)*a3*k*s1*(s2*s2)*s3*1.6E1-(a*a*a*a)*a3*k*(s1*s1)*s2*s3*4.0-a*a2*la*s1*(s2*s2)*(s3*s3)*2.0-a*a2*la*(s1*s1)*s2*(s3*s3)*2.0-a*a2*la*(s1*s1)*(s2*s2)*s3-(a*a)*a2*la*s1*s2*(s3*s3)*1.2E1-(a*a)*a2*la*s1*(s2*s2)*s3*4.0-(a*a)*a2*la*(s1*s1)*s2*s3*4.0+a*a3*la*s1*(s2*s2)*(s3*s3)*2.0+a*a3*la*(s1*s1)*s2*(s3*s3)+a*a3*la*(s1*s1)*(s2*s2)*s3*2.0+(a*a)*a3*la*s1*s2*(s3*s3)*4.0+(a*a)*a3*la*s1*(s2*s2)*s3*1.2E1+(a*a)*a3*la*(s1*s1)*s2*s3*4.0-a*a2*mu*(s1*s1)*(s2*s2)*s3-(a*a)*a2*mu*s1*(s2*s2)*s3*2.0-(a*a)*a2*mu*(s1*s1)*s2*s3*2.0-(a*a*a)*a2*mu*s1*s2*(s3*s3*s3)*1.2E1-(a*a*a*a)*a2*mu*s1*s2*(s3*s3)*2.4E1-(a*a*a*a)*a2*mu*s1*(s2*s2)*s3*6.0-(a*a*a*a)*a2*mu*(s1*s1)*s2*s3*6.0+a*a3*mu*(s1*s1)*s2*(s3*s3)+(a*a)*a3*mu*s1*s2*(s3*s3)*2.0+(a*a)*a3*mu*(s1*s1)*s2*s3*2.0+(a*a*a)*a3*mu*s1*(s2*s2*s2)*s3*1.2E1+(a*a*a*a)*a3*mu*s1*s2*(s3*s3)*6.0+(a*a*a*a)*a3*mu*s1*(s2*s2)*s3*2.4E1+(a*a*a*a)*a3*mu*(s1*s1)*s2*s3*6.0+a*la*(s1*s1)*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-a*la*(s1*s1)*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*la*s1*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a)*la*s1*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*4.0+(a*a)*a2*la*(s1*s1)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a3*la*(s1*s1)*(s2*s2)*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*a2*la*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a3*la*(s2*s2)*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0+a*a2*la*(s1*s1)*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))+(a*a)*a2*la*s1*(s2*s2)*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a)*a2*la*(s1*s1)*s2*s3*log((a+s1)*(a+s2)*(a+s3))*2.0-a*a3*la*(s1*s1)*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))-(a*a)*a3*la*s1*s2*(s3*s3)*log((a+s1)*(a+s2)*(a+s3))*2.0-(a*a)*a3*la*(s1*s1)*s2*s3*log((a+s1)*(a+s2)*(a+s3))*2.0+(a*a*a)*a2*la*s1*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0-(a*a*a)*a3*la*s1*s2*s3*log((a+s1)*(a+s2)*(a+s3))*4.0))/(s2*s2-s3*s3);
    
    
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
