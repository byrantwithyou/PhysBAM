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
    
    a1 = (a*a-a*(s2+s3)+s2*s3)/den;
    a2 = (a*a-a*(s1+s3)+s1*s3)/den;
    a3 = (a*a-a*(s2+s1)+s2*s1)/den;
    a11 = -2.0*a1*(s3-2.0*a-a1*s3+3.0*a1*a+s2-a1*s2-a1*s1)/den;
    a22 = -2.0*a2*(s3-2.0*a-a2*s3+3.0*a2*a+s1-a2*s1-a2*s2)/den;
    a33 = -2.0*a3*(s1-2.0*a-a3*s1+3.0*a3*a+s2-a3*s2-a3*s3)/den;
    a12 = -(-s3+a+a1*s3-2.0*a1*a+a1*s1-2.0*a1*a2*s3+6.0*a1*a2*a-2.0*a1*a2*s2-2.0*a1*a2*s1+a2*s3-2.0*a2*a+a2*s2)/den;
    a13 = -(-s2+a+a1*s2-2.0*a1*a+a1*s1-2.0*a1*a3*s2+6.0*a1*a3*a-2.0*a1*a3*s3-2.0*a1*a3*s1+a3*s2-2.0*a3*a+a3*s3)/den;
    a23 = -(-s1+a+a3*s1-2.0*a3*a+a3*s3-2.0*a3*a2*s1+6.0*a3*a2*a-2.0*a3*a2*s2-2.0*a3*a2*s3+a2*s1-2.0*a2*a+a2*s2)/den;
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
    
    a1 = (a*a-a*(s2+s3)+s2*s3)/den;
    a2 = (a*a-a*(s1+s3)+s1*s3)/den;
    a3 = (a*a-a*(s2+s1)+s2*s1)/den;
    a11 = -2.0*a1*(s3-2.0*a-a1*s3+3.0*a1*a+s2-a1*s2-a1*s1)/den;
    a22 = -2.0*a2*(s3-2.0*a-a2*s3+3.0*a2*a+s1-a2*s1-a2*s2)/den;
    a33 = -2.0*a3*(s1-2.0*a-a3*s1+3.0*a3*a+s2-a3*s2-a3*s3)/den;
    a12 = -(-s3+a+a1*s3-2.0*a1*a+a1*s1-2.0*a1*a2*s3+6.0*a1*a2*a-2.0*a1*a2*s2-2.0*a1*a2*s1+a2*s3-2.0*a2*a+a2*s2)/den;
    a13 = -(-s2+a+a1*s2-2.0*a1*a+a1*s1-2.0*a1*a3*s2+6.0*a1*a3*a-2.0*a1*a3*s3-2.0*a1*a3*s1+a3*s2-2.0*a3*a+a3*s3)/den;
    a23 = -(-s1+a+a3*s1-2.0*a3*a+a3*s3-2.0*a3*a2*s1+6.0*a3*a2*a-2.0*a3*a2*s2-2.0*a3*a2*s3+a2*s1-2.0*a2*a+a2*s2)/den;
    /////////////Hopefully I can get the member function to work at some point    
    
    
    T t1,t2,t3,t4,t5,t6,t7,t13,t15,t19,t21,t26,t32,t37;

    T mu = constant_mu;
    T la = constant_lambda;
    //Calculate_A(s1,s2,s3);
    T k = extra_force_coefficient; 
    
    t1 = a;
    t2 = s1 + t1;
    t3 = (int) (t2 * t2);
    t4 = s2 + t1;
    t5 = (int) (t4 * t4);
    t6 = s3 + t1;
    t7 = (int) (t6 * t6);
    t13 = log(t6 * t4 * t2);
    t15 = t13 * t13;
    t19 = 0.1e1 / t2;
    t21 = t13 * la;
    t26 = 0.1e1 / t4;
    t32 = 0.1e1 / t6;
    t37 = (double) ((t3 + t5 + t7 - 2) * mu) / 0.2e1 - (double) mu * t13 + la * t15 / 0.2e1 - (t2 * (double) mu - (double) mu * t19 + t19 * t21) * t1 - (t4 * (double) mu - t26 * (double) mu + t26 * t21) * t1 - (t6 * (double) mu - t32 * (double) mu + t32 * t21) * t1 +k*a*a;

    return t37;
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
    
    a1 = (a*a-a*(s2+s3)+s2*s3)/den;
    a2 = (a*a-a*(s1+s3)+s1*s3)/den;
    a3 = (a*a-a*(s2+s1)+s2*s1)/den;
    a11 = -2.0*a1*(s3-2.0*a-a1*s3+3.0*a1*a+s2-a1*s2-a1*s1)/den;
    a22 = -2.0*a2*(s3-2.0*a-a2*s3+3.0*a2*a+s1-a2*s1-a2*s2)/den;
    a33 = -2.0*a3*(s1-2.0*a-a3*s1+3.0*a3*a+s2-a3*s2-a3*s3)/den;
    a12 = -(-s3+a+a1*s3-2.0*a1*a+a1*s1-2.0*a1*a2*s3+6.0*a1*a2*a-2.0*a1*a2*s2-2.0*a1*a2*s1+a2*s3-2.0*a2*a+a2*s2)/den;
    a13 = -(-s2+a+a1*s2-2.0*a1*a+a1*s1-2.0*a1*a3*s2+6.0*a1*a3*a-2.0*a1*a3*s3-2.0*a1*a3*s1+a3*s2-2.0*a3*a+a3*s3)/den;
    a23 = -(-s1+a+a3*s1-2.0*a3*a+a3*s3-2.0*a3*a2*s1+6.0*a3*a2*a-2.0*a3*a2*s2-2.0*a3*a2*s3+a2*s1-2.0*a2*a+a2*s2)/den;
    /////////////Hopefully I can get the member function to work at some point
    
    T mu = constant_mu;
    T la = constant_lambda;
    DIAGONAL_MATRIX<T,3> result;
    //Calculate_A(s1,s2,s3);
    T k = extra_force_coefficient;    
    T t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42;
    
    
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
    result.x33 = -a3*t21-a3*t24+mu*(a3*t26+a3*t28+t30*t37)*(1.0/2.0)-a*(t38+a3*mu*t12-a3*la*t12*t8+la*t12*t13*t42*t9)-a*(t38+a3*mu*t4-a3*la*t4*t8+la*t10*t13*t4*t42)-mu*t10*t13*t42*t9+la*t10*t13*t42*t8*t9+2.0*k*a3*a;
    
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
    
    a1 = (a*a-a*(s2+s3)+s2*s3)/den;
    a2 = (a*a-a*(s1+s3)+s1*s3)/den;
    a3 = (a*a-a*(s2+s1)+s2*s1)/den;
    a11 = -2.0*a1*(s3-2.0*a-a1*s3+3.0*a1*a+s2-a1*s2-a1*s1)/den;
    a22 = -2.0*a2*(s3-2.0*a-a2*s3+3.0*a2*a+s1-a2*s1-a2*s2)/den;
    a33 = -2.0*a3*(s1-2.0*a-a3*s1+3.0*a3*a+s2-a3*s2-a3*s3)/den;
    a12 = -(-s3+a+a1*s3-2.0*a1*a+a1*s1-2.0*a1*a2*s3+6.0*a1*a2*a-2.0*a1*a2*s2-2.0*a1*a2*s1+a2*s3-2.0*a2*a+a2*s2)/den;
    a13 = -(-s2+a+a1*s2-2.0*a1*a+a1*s1-2.0*a1*a3*s2+6.0*a1*a3*a-2.0*a1*a3*s3-2.0*a1*a3*s1+a3*s2-2.0*a3*a+a3*s3)/den;
    a23 = -(-s1+a+a3*s1-2.0*a3*a+a3*s3-2.0*a3*a2*s1+6.0*a3*a2*a-2.0*a3*a2*s2-2.0*a3*a2*s3+a2*s1-2.0*a2*a+a2*s2)/den;
    /////////////Hopefully I can get the member function to work at some point
    
    T t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87,t88,t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150,t151,t152,t153,t154,t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187,t188,t189,t190,t191,t192,t193,t194,t195,t196,t197,t198,t199,t200,t201,t202,t203,t204,t205,t206,t207,t208,t209,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t224,t225,t226,t227,t228,t229,t230,t231,t232,t233,t234,t235,t236,t237,t238,t239,t240,t241,t242,t243,t244,t245,t246,t247,t248,t249,t250,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262;
    T b1,b2,b3;
    
    t44 = a1+1.0;
    t45 = a+s1;
    t46 = 1.0/(t45*t45);
    t47 = a+s2;
    t48 = a+s3;
    t49 = t45*t47*t48;
    t50 = log(t49);
    t51 = 1.0/t45;
    t52 = 1.0/t47;
    t53 = a*2.0;
    t54 = a1*a1;
    t55 = 1.0/(t47*t47);
    t56 = 1.0/(t47*t47*t47);
    t57 = 1.0/t48;
    t58 = s1*2.0;
    t59 = t53+t58;
    t60 = s2*2.0;
    t61 = t53+t60;
    t62 = s3*2.0;
    t63 = t53+t62;
    t64 = a1*t45*t47;
    t65 = a1*t45*t48;
    t66 = t44*t47*t48;
    t67 = t64+t65+t66;
    t68 = t44*t44;
    t69 = 1.0/(t48*t48);
    t70 = 1.0/(t45*t45*t45);
    t71 = t54*t59;
    t72 = a11*t45*t47;
    t73 = a11*t45*t48;
    t74 = a11*t47*t48;
    t75 = a1*t44*t61;
    t76 = a1*t44*t63;
    t77 = t71+t72+t73+t74+t75+t76;
    t78 = a2+1.0;
    t79 = mu*t45;
    t80 = la*t50*t51;
    t98 = mu*t51;
    t81 = t79+t80-t98;
    t82 = mu*t47;
    t83 = la*t50*t52;
    t99 = mu*t52;
    t84 = t82+t83-t99;
    t85 = a2*a2;
    t86 = a2*t45*t47;
    t87 = a2*t47*t48;
    t88 = t45*t48*t78;
    t89 = t86+t87+t88;
    t90 = t78*t78;
    t91 = t61*t85;
    t92 = a22*t45*t47;
    t93 = a22*t45*t48;
    t94 = a22*t47*t48;
    t95 = a2*t59*t78;
    t96 = a2*t63*t78;
    t97 = t91+t92+t93+t94+t95+t96;
    t100 = a3+1.0;
    t101 = a3*a3;
    t102 = a3*t45*t48;
    t103 = a3*t47*t48;
    t104 = t100*t45*t47;
    t105 = t102+t103+t104;
    t106 = a3*mu;
    t107 = t101*t63;
    t108 = a33*t45*t47;
    t109 = a33*t45*t48;
    t110 = a33*t47*t48;
    t111 = a3*t100*t59;
    t112 = a3*t100*t61;
    t113 = t107+t108+t109+t110+t111+t112;
    t114 = mu*t44;
    t115 = mu*t44*t46;
    t116 = la*t46*t52*t57*t67;
    t140 = la*t44*t46*t50;
    t117 = t114+t115+t116-t140;
    t118 = mu*t78;
    t119 = mu*t55*t78;
    t120 = la*t51*t55*t57*t89;
    t160 = la*t50*t55*t78;
    t121 = t118+t119+t120-t160;
    t122 = a1*mu;
    t123 = a1*mu*t55;
    t124 = la*t51*t55*t57*t67;
    t143 = a1*la*t50*t55;
    t125 = t122+t123+t124-t143;
    t126 = a2*mu;
    t127 = a2*mu*t46;
    t128 = la*t46*t52*t57*t89;
    t165 = a2*la*t46*t50;
    t129 = t126+t127+t128-t165;
    t130 = t44*t48*t78;
    t131 = a12*t45*t47;
    t132 = a12*t45*t48;
    t133 = a12*t47*t48;
    t134 = a1*a2*t45;
    t135 = a1*a2*t47;
    t136 = a1*a2*t48;
    t137 = a1*t45*t78;
    t138 = a2*t44*t47;
    t139 = t130+t131+t132+t133+t134+t135+t136+t137+t138;
    t141 = a1*2.0;
    t142 = t141+2.0;
    t144 = a3*mu*t46;
    t145 = la*t105*t46*t52*t57;
    t166 = a3*la*t46*t50;
    t146 = t106+t144+t145-t166;
    t147 = a3*mu*t55;
    t148 = la*t105*t51*t55*t57;
    t167 = a3*la*t50*t55;
    t149 = t106+t147+t148-t167;
    t150 = t100*t44*t47;
    t151 = a13*t45*t47;
    t152 = a13*t45*t48;
    t153 = a13*t47*t48;
    t154 = a1*a3*t45;
    t155 = a1*a3*t47;
    t156 = a1*a3*t48;
    t157 = a1*t100*t45;
    t158 = a3*t44*t48;
    t159 = t150+t151+t152+t153+t154+t155+t156+t157+t158;
    t161 = a3*2.0;
    t162 = t161+2.0;
    t163 = a2*2.0;
    t164 = t163+2.0;
    t168 = t100*t45*t78;
    t169 = a23*t45*t47;
    t170 = a23*t45*t48;
    t171 = a23*t47*t48;
    t172 = a2*a3*t45;
    t173 = a2*a3*t47;
    t174 = a2*a3*t48;
    t175 = a2*t100*t47;
    t176 = a3*t48*t78;
    t177 = t168+t169+t170+t171+t172+t173+t174+t175+t176;
    t178 = a*a;
    t179 = t178*t178;
    t180 = s1*s1;
    t181 = s2*s2;
    t182 = t181*t181;
    t183 = t180*t180;
    t184 = s3*s3;
    t185 = a*mu*t178*t180;
    t186 = a*mu*s1*s3*t178*2.0;
    t187 = s1*t178;
    t188 = s2*t178;
    t189 = s3*t178;
    t190 = a*s1*s2;
    t191 = a*s1*s3;
    t192 = a*s2*s3;
    t193 = s1*s2*s3;
    t194 = t187+t188+t189+t190+t191+t192+t193+a*t178;
    t195 = log(t194);
    t196 = mu*s3*t179;
    t197 = mu*s3*t178*t184;
    t198 = a*mu*t178*t184*2.0;
    t199 = a3*mu*s3*t178*t184;
    t200 = la*s3*t178*t195;
    t201 = a*a3*la*t184*t195;
    t202 = a3*la*s3*t178*t195*2.0;
    t240 = la*t195;
    t203 = la+mu-t240;
    t204 = t180-t181; if(fabs(t204)<panic_threshold) t204=t204<0?-panic_threshold:panic_threshold;
    t205 = 1.0/t204;
    t206 = mu*s1*t179*2.0;
    t207 = la*s1*t179*2.0;
    t208 = la*s2*t179*t50*2.0;
    t209 = a*la*s1*s3*t178*2.0;
    t210 = a*mu*s2*s3*t178*2.0;
    t211 = a*a1*la*t180*t181;
    t212 = a*a1*mu*t180*t181*2.0;
    t213 = a*a1*mu*t178*t180*t181*2.0;
    t214 = a*a1*la*s1*s2*t178*7.0;
    t215 = a*a1*mu*s1*s2*t178*4.0;
    t216 = a*a1*mu*s1*s2*t179*2.0;
    t217 = a*la*s2*s3*t178*t50*2.0;
    t218 = a*a2*mu*t180*t181*t184*2.0;
    t219 = a*a2*la*s1*s2*t178*t50*4.0;
    t220 = a1*la*s1*s2*s3*t178*4.0;
    t221 = a1*mu*s1*s2*s3*t178*2.0;
    t222 = a*a2*la*t180*t181*t50*2.0;
    t223 = a*a2*mu*s1*s2*t178*t184*2.0;
    t224 = a2*la*s1*s2*s3*t178*t50*2.0;
    t225 = a*t183;
    t226 = s3*t183;
    t227 = s1*t178*t180*2.0;
    t228 = a*t178*t180;
    t229 = a*s1*s3*t180*2.0;
    t230 = s3*t178*t180;
    t247 = s3*t178*t184;
    t248 = a*t178*t184;
    t231 = t225+t226+t227+t228+t229+t230-t247-t248-a*t180*t184-s1*t178*t184*2.0-s3*t180*t184-a*s1*s3*t184*2.0; if(fabs(t231)<panic_threshold) t231=t231<0?-panic_threshold:panic_threshold;
    t232 = 1.0/t231;
    t233 = la*s1*t178;
    t234 = a*a1*la*s1*s3;
    t235 = a*a1*mu*s1*s3*t178*2.0;
    t236 = s1+s3;
    t237 = 1.0/t236;if(fabs(t237)<panic_threshold) t237=t237<0?-panic_threshold:panic_threshold;
    t238 = a*s1;
    t239 = t45*t45;
    t241 = s1-s3;if(fabs(t241)<panic_threshold) t241=t241<0?-panic_threshold:panic_threshold;
    t242 = 1.0/t241;
    t243 = a*t182;
    t244 = s3*t182;
    t245 = s2*t178*t181*2.0;
    t246 = a*t178*t181;
    t249 = a*s2*s3*t181*2.0;
    t250 = s3*t178*t181;
    t251 = a*mu*t184*2.0;
    t252 = mu*s3*t178*2.0;
    t253 = la*s2*t178;
    t254 = a*la*t184;
    t255 = la*s3*t178;
    t256 = a*a2*la*s2*s3;
    t257 = s2+s3;if(fabs(t257)<panic_threshold) t257=t257<0?-panic_threshold:panic_threshold;
    t258 = 1.0/t257;
    t259 = a*s2;
    t260 = t47*t47;
    t261 = s2-s3;if(fabs(t261)<panic_threshold) t261=t261<0?-panic_threshold:panic_threshold;
    t262 = 1.0/t261;
    b1 = s1*s1-s2*s2;if(fabs(b1)<panic_threshold) b1=b1<0?-panic_threshold:panic_threshold;
    b2 = s1*s1-s3*s3;if(fabs(b2)<panic_threshold) b2=b2<0?-panic_threshold:panic_threshold;
    b3 = s2*s2-s3*s3;if(fabs(b3)<panic_threshold) b3=b3<0?-panic_threshold:panic_threshold;
    dP_dF.x1111 = a1*t117*-2.0-a1*t125*2.0-a11*t81-a11*t84+mu*(t54*4.0+t68*2.0+a11*t59+a11*t61+a11*t63)*(1.0/2.0)+a*(-a11*mu-a11*mu*t55+mu*t54*t56*2.0+a11*la*t50*t55-la*t50*t54*t56*2.0-la*t51*t55*t57*t77+la*t44*t46*t55*t57*t67+a1*la*t51*t56*t57*t67*3.0+a1*la*t51*t55*t67*t69)-a*a11*mu-a*a11*mu*t46+a*mu*t68*t70*2.0+a*a11*la*t46*t50-a*la*t50*t68*t70*2.0-mu*t51*t52*t57*t77+la*t46*t55*(t67*t67)*t69+la*t50*t51*t52*t57*t77+mu*t44*t46*t52*t57*t67-a*la*t46*t52*t57*t77+a1*mu*t51*t55*t57*t67+a1*mu*t51*t52*t67*t69+a*a1*la*t46*t55*t57*t67+a*a1*la*t46*t52*t67*t69+a*la*t44*t52*t57*t67*t70*3.0-a1*la*t50*t51*t55*t57*t67-a1*la*t50*t51*t52*t67*t69-la*t44*t46*t50*t52*t57*t67+2.0*k*a*a11 + 2.0*k*a1*a1;
    dP_dF.x2222 = a2*t121*-2.0-a2*t129*2.0-a22*t81-a22*t84+mu*(t85*4.0+t90*2.0+a22*t59+a22*t61+a22*t63)*(1.0/2.0)+a*(-a22*mu-a22*mu*t46+mu*t70*t85*2.0+a22*la*t46*t50-la*t50*t70*t85*2.0-la*t46*t52*t57*t97+la*t46*t55*t57*t78*t89+a2*la*t46*t52*t69*t89+a2*la*t52*t57*t70*t89*3.0)-a*a22*mu-a*a22*mu*t55+a*mu*t56*t90*2.0+a*a22*la*t50*t55-a*la*t50*t56*t90*2.0-mu*t51*t52*t57*t97+la*t46*t55*t69*(t89*t89)+la*t50*t51*t52*t57*t97+mu*t51*t55*t57*t78*t89-a*la*t51*t55*t57*t97+a2*mu*t46*t52*t57*t89+a2*mu*t51*t52*t69*t89+a*a2*la*t46*t55*t57*t89+a*a2*la*t51*t55*t69*t89+a*la*t51*t56*t57*t78*t89*3.0-a2*la*t46*t50*t52*t57*t89-a2*la*t50*t51*t52*t69*t89-la*t50*t51*t55*t57*t78*t89+2.0*k*a*a22 + 2.0*k*a2*a2;
    dP_dF.x3333 = a3*t146*-2.0-a3*t149*2.0-a33*t81-a33*t84+a*(-a33*mu-a33*mu*t55+mu*t101*t56*2.0+a33*la*t50*t55-la*t101*t50*t56*2.0-la*t113*t51*t55*t57+la*t100*t105*t51*t55*t69+a3*la*t105*t46*t55*t57+a3*la*t105*t51*t56*t57*3.0)+mu*(t101*4.0+a33*t59+a33*t61+a33*t63+(t100*t100)*2.0)*(1.0/2.0)-a*a33*mu-a*a33*mu*t46+a*mu*t101*t70*2.0+a*a33*la*t46*t50-a*la*t101*t50*t70*2.0-mu*t113*t51*t52*t57+la*(t105*t105)*t46*t55*t69+la*t113*t50*t51*t52*t57+mu*t100*t105*t51*t52*t69-a*la*t113*t46*t52*t57+a3*mu*t105*t46*t52*t57+a3*mu*t105*t51*t55*t57+a*a3*la*t105*t46*t55*t57+a*a3*la*t105*t52*t57*t70*3.0+a*la*t100*t105*t46*t52*t69-a3*la*t105*t46*t50*t52*t57-a3*la*t105*t50*t51*t55*t57-la*t100*t105*t50*t51*t52*t69+2.0*k*a*a33 + 2.0*k*a3*a3;
    dP_dF.x2211 = -a1*t121-a1*t129-a2*t117-a2*t125-a12*t81-a12*t84+a*(-a12*mu-a12*mu*t55+a12*la*t50*t55+a1*mu*t56*t78*2.0-a1*la*t50*t56*t78*2.0-la*t139*t51*t55*t57+la*t51*t56*t57*t67*t78*2.0+a2*la*t46*t55*t57*t67+a2*la*t51*t55*t67*t69+a1*la*t51*t56*t57*t89)+mu*(a1*a2*2.0+a1*t164+a2*t142+a12*t59+a12*t61+a12*t63)*(1.0/2.0)-a*a12*mu-a*a12*mu*t46+a*a12*la*t46*t50+a*a2*mu*t44*t70*2.0-mu*t139*t51*t52*t57+la*t139*t50*t51*t52*t57+la*t46*t55*t67*t69*t89+mu*t51*t55*t57*t67*t78-a*a2*la*t44*t50*t70*2.0-a*la*t139*t46*t52*t57+a2*mu*t46*t52*t57*t67+a2*mu*t51*t52*t67*t69+a*a2*la*t46*t52*t67*t69+a*a2*la*t52*t57*t67*t70*2.0+a*la*t46*t55*t57*t67*t78+a*la*t44*t52*t57*t70*t89-a2*la*t46*t50*t52*t57*t67-a2*la*t50*t51*t52*t67*t69-la*t50*t51*t55*t57*t67*t78+2.0*k*a*a12 + 2.0*k*a2*a1;
    dP_dF.x3311 = a*(-a13*mu-a13*mu*t55+a1*a3*mu*t56*2.0+a13*la*t50*t55-a1*a3*la*t50*t56*2.0-la*t159*t51*t55*t57+la*t100*t51*t55*t67*t69+a1*la*t105*t51*t56*t57+a3*la*t46*t55*t57*t67+a3*la*t51*t56*t57*t67*2.0)-a1*t146-a1*t149-a3*t117-a3*t125-a13*t81-a13*t84+mu*(a1*a3*2.0+a1*t162+a3*t142+a13*t59+a13*t61+a13*t63)*(1.0/2.0)-a*a13*mu-a*a13*mu*t46+a*a13*la*t46*t50+a*a3*mu*t44*t70*2.0-mu*t159*t51*t52*t57+la*t159*t50*t51*t52*t57+la*t105*t46*t55*t67*t69+mu*t100*t51*t52*t67*t69-a*a3*la*t44*t50*t70*2.0-a*la*t159*t46*t52*t57+a3*mu*t46*t52*t57*t67+a3*mu*t51*t55*t57*t67+a*a3*la*t46*t55*t57*t67+a*a3*la*t52*t57*t67*t70*2.0+a*la*t105*t44*t52*t57*t70+a*la*t100*t46*t52*t67*t69-a3*la*t46*t50*t52*t57*t67-a3*la*t50*t51*t55*t57*t67-la*t100*t50*t51*t52*t67*t69+2.0*k*a*a13 + 2.0*k*a3*a1;
    dP_dF.x3322 = a*(-a23*mu-a23*mu*t46+a2*a3*mu*t70*2.0+a23*la*t46*t50-a2*a3*la*t50*t70*2.0-la*t177*t46*t52*t57+la*t100*t46*t52*t69*t89+a2*la*t105*t52*t57*t70+a3*la*t46*t55*t57*t89+a3*la*t52*t57*t70*t89*2.0)-a2*t146-a2*t149-a3*t121-a3*t129-a23*t81-a23*t84+mu*(a2*a3*2.0+a2*t162+a3*t164+a23*t59+a23*t61+a23*t63)*(1.0/2.0)-a*a23*mu-a*a23*mu*t55+a*a23*la*t50*t55+a*a3*mu*t56*t78*2.0-mu*t177*t51*t52*t57+la*t177*t50*t51*t52*t57+la*t105*t46*t55*t69*t89+mu*t100*t51*t52*t69*t89-a*a3*la*t50*t56*t78*2.0-a*la*t177*t51*t55*t57+a3*mu*t46*t52*t57*t89+a3*mu*t51*t55*t57*t89+a*a3*la*t46*t55*t57*t89+a*a3*la*t51*t56*t57*t89*2.0+a*la*t105*t51*t56*t57*t78+a*la*t100*t51*t55*t69*t89-a3*la*t46*t50*t52*t57*t89-a3*la*t50*t51*t55*t57*t89-la*t100*t50*t51*t52*t69*t89+2.0*k*a*a23 + 2.0*k*a2*a3;
    dP_dF.x2121 = -t205*t46*t55*t57*(t185+t186+t206+t207+t208+t209+t211+t212+t213+t214+t215+t216+t217+t218+t219+t220+t221+t222+t223+t224+t235-la*s2*t179*2.0-mu*s2*t179*2.0+a1*la*s1*t179*6.0-a2*la*s2*t179*6.0+a*la*t178*t180-a*la*t178*t181+a1*mu*s1*t179*3.0-a2*mu*s2*t179*3.0-a*mu*t178*t181-a*mu*t179*t180+a*mu*t178*t182+a*mu*t179*t181-a*mu*t178*t183+a*mu*t180*t182-a*mu*t181*t183+la*s3*t178*t180-la*s3*t178*t181-la*s1*t179*t50*2.0-mu*s1*t179*t180*2.0+mu*s1*t178*t182*2.0+mu*s1*t179*t181*2.0-mu*s2*t179*t180*2.0+mu*s2*t179*t181*2.0-mu*s2*t178*t183*2.0+mu*s3*t178*t180-mu*s3*t178*t181-mu*s3*t179*t180+mu*s3*t178*t182+mu*s3*t179*t181-mu*s3*t178*t183+mu*s3*t180*t182-mu*s3*t181*t183+a*a1*la*t178*t180*7.0-a*a2*la*t178*t181*7.0-a*a2*la*t180*t181+a*a1*mu*t178*t180*4.0+a*a1*mu*t179*t180*2.0-a*a2*mu*t178*t181*4.0-a*a2*mu*t179*t181*2.0-a*a2*mu*t180*t181*2.0-a*la*s2*s3*t178*2.0+a*mu*s1*s3*t182*2.0-a*mu*s2*s3*t178*2.0-a*mu*s2*s3*t183*2.0+a1*la*s1*t178*t180*2.0+a1*la*s1*t178*t181*2.0+a1*la*s2*t178*t180*6.0-a2*la*s1*t178*t181*6.0+a1*la*s3*t178*t180*4.0-a2*la*s2*t178*t180*2.0-a2*la*s2*t178*t181*2.0-a2*la*s3*t178*t181*4.0-a1*la*s1*t179*t50*3.0+a2*la*s2*t179*t50*3.0-a*la*t178*t180*t50+a*la*t178*t181*t50+a1*mu*s1*t178*t179+a1*mu*s1*t178*t180*2.0+a1*mu*s1*t178*t181*2.0+a1*mu*s1*t179*t180+a1*mu*s1*t179*t181-a1*mu*s1*t179*t184+a1*mu*s1*t180*t181+a1*mu*s2*t178*t180*4.0+a1*mu*s2*t179*t180*4.0-a2*mu*s1*t178*t181*4.0-a2*mu*s1*t179*t181*4.0-a2*mu*s2*t178*t179+a1*mu*s3*t178*t180*2.0-a2*mu*s2*t178*t180*2.0-a2*mu*s2*t178*t181*2.0-a2*mu*s2*t179*t180-a2*mu*s2*t179*t181+a2*mu*s2*t179*t184-a2*mu*s2*t180*t181-a2*mu*s3*t178*t181*2.0-la*s3*t178*t180*t50+la*s3*t178*t181*t50-mu*s1*t178*t180*t181*2.0+mu*s2*t178*t180*t181*2.0-mu*s1*s2*s3*t178*t180*4.0+mu*s1*s2*s3*t178*t181*4.0+a*a1*la*s1*s2*t180+a*a1*la*s1*s3*t178*4.0-a*a2*la*s1*s2*t178*7.0+a*a1*la*s1*s3*t180+a*a1*la*s1*s3*t181-a*a2*la*s1*s2*t181+a*a1*la*s2*s3*t180*2.0-a*a2*la*s1*s3*t181*2.0-a*a2*la*s2*s3*t178*4.0-a*a2*la*s2*s3*t180-a*a2*la*s2*s3*t181+a*a1*mu*s1*s2*t180*2.0-a*a2*mu*s1*s2*t178*4.0-a*a2*mu*s1*s2*t179*2.0+a*a1*mu*s1*s3*t180+a*a1*mu*s1*s3*t181-a*a2*mu*s1*s2*t181*2.0-a*a2*mu*s2*s3*t178*2.0-a*a2*mu*s2*s3*t180-a*a2*mu*s2*s3*t181-a*a1*la*t178*t180*t50*4.0-a*a1*la*t180*t181*t50*2.0+a*a2*la*t178*t181*t50*4.0-a*a1*mu*t178*t180*t184*2.0-a*a1*mu*t180*t181*t184*2.0-a*a2*mu*t178*t180*t181*2.0+a*a2*mu*t178*t181*t184*2.0-a2*la*s1*s2*s3*t178*4.0-a*la*s1*s3*t178*t50*2.0-a2*mu*s1*s2*s3*t178*2.0-a*mu*s1*s2*t178*t180*4.0+a*mu*s1*s2*t178*t181*4.0-a*mu*s1*s3*t178*t180*2.0+a*mu*s1*s3*t178*t181*2.0-a*mu*s1*s3*t180*t181*2.0-a*mu*s2*s3*t178*t180*2.0+a*mu*s2*s3*t178*t181*2.0+a*mu*s2*s3*t180*t181*2.0-a1*la*s1*t178*t180*t50*2.0-a1*la*s1*t178*t181*t50*2.0-a1*la*s1*t180*t181*t50-a1*la*s2*t178*t180*t50*4.0+a2*la*s1*t178*t181*t50*4.0-a1*la*s3*t178*t180*t50*2.0+a2*la*s2*t178*t180*t50*2.0+a2*la*s2*t178*t181*t50*2.0+a2*la*s2*t180*t181*t50+a2*la*s3*t178*t181*t50*2.0+a1*mu*s1*t178*t180*t181-a1*mu*s1*t178*t180*t184-a1*mu*s1*t178*t181*t184-a1*mu*s1*t180*t181*t184-a1*mu*s2*t178*t180*t184*4.0+a2*mu*s1*t178*t181*t184*4.0-a2*mu*s2*t178*t180*t181+a2*mu*s2*t178*t180*t184+a2*mu*s2*t178*t181*t184+a2*mu*s2*t180*t181*t184-a*a1*la*s1*s2*t178*t50*4.0-a*a1*la*s1*s2*t180*t50*2.0-a*a1*la*s1*s3*t178*t50*2.0-a*a1*la*s1*s3*t180*t50-a*a1*la*s1*s3*t181*t50+a*a2*la*s1*s2*t181*t50*2.0+a*a2*la*s2*s3*t178*t50*2.0+a*a2*la*s2*s3*t180*t50+a*a2*la*s2*s3*t181*t50+a*a1*mu*s1*s2*t178*t180*2.0-a*a1*mu*s1*s2*t178*t184*2.0-a*a1*mu*s1*s2*t180*t184*2.0-a*a2*mu*s1*s2*t178*t181*2.0+a*a2*mu*s1*s2*t181*t184*2.0-a1*la*s1*s2*s3*t178*t50*2.0)+(2.0*k*(s1*a1 - s2*a2)*a)/(b1);
    dP_dF.x3131 = -t232*(-t185+t186+t196+t197+t198+t199+t200+t201+t202+t233+t234+a*mu*t180-a*mu*t183-la*s3*t178+mu*s1*t178*2.0-mu*s3*t178-mu*s3*t183+a*a1*la*t180-a*a3*la*t184+a*a1*mu*t180*2.0-a*a3*mu*t184+a1*la*s1*t178*2.0-a3*la*s3*t178*2.0-a*la*t180*t195+a1*mu*s1*t178*2.0+a1*mu*s1*t179+a1*mu*s1*t180-a3*mu*s3*t178*2.0-a3*mu*s3*t179-a3*mu*s3*t180+a*mu*t180*t184*2.0-la*s1*t178*t195*2.0-mu*s1*t178*t180*2.0+mu*s1*t178*t184*4.0+mu*s3*t180*t184-a*a3*la*s1*s3+a*a1*mu*s1*s3-a*a3*mu*s1*s3*2.0-a*a1*la*t180*t195*2.0+a*a1*mu*t178*t180*2.0-a*a1*mu*t180*t184*2.0-a*mu*s1*s3*t180*2.0+a*mu*s1*s3*t184*2.0-a1*la*s1*t178*t195*2.0-a1*la*s1*t180*t195+a3*la*s3*t180*t195+a1*mu*s1*t178*t180-a1*mu*s1*t178*t184-a1*mu*s1*t180*t184-a3*mu*s3*t178*t180+a3*mu*s3*t180*t184-a*a1*la*s1*s3*t195+a*a3*la*s1*s3*t195*2.0-a*a3*mu*s1*s3*t178*2.0+a*a3*mu*s1*s3*t184*2.0)-t237*t242*t46*t57*(a*la*t45*t52*(t238-a*s3+a1*t180-a3*t184*2.0+a*a1*s1*3.0-a*a3*s3*3.0+a1*s1*s3*2.0-a3*s1*s3)+a*t203*t239*t48*t55*(a1*s1-a3*s3))+(2.0*k*(s1*a1 - s3*a3)*a)/(b2);
    dP_dF.x3232 = -(t196+t197+t198+t199+t200+t201+t202+t210+t253+t256+a*mu*t181-a*mu*t182-la*s3*t178+mu*s2*t178*2.0-mu*s3*t178-mu*s3*t182+a*a2*la*t181-a*a3*la*t184+a*a2*mu*t181*2.0-a*a3*mu*t184+a2*la*s2*t178*2.0-a3*la*s3*t178*2.0-a*la*t181*t195+a2*mu*s2*t178*2.0+a2*mu*s2*t179+a2*mu*s2*t181-a3*mu*s3*t178*2.0-a3*mu*s3*t179-a3*mu*s3*t181-a*mu*t178*t181+a*mu*t181*t184*2.0-la*s2*t178*t195*2.0-mu*s2*t178*t181*2.0+mu*s2*t178*t184*4.0+mu*s3*t181*t184-a*a3*la*s2*s3+a*a2*mu*s2*s3-a*a3*mu*s2*s3*2.0-a*a2*la*t181*t195*2.0+a*a2*mu*t178*t181*2.0-a*a2*mu*t181*t184*2.0-a*mu*s2*s3*t181*2.0+a*mu*s2*s3*t184*2.0-a2*la*s2*t178*t195*2.0-a2*la*s2*t181*t195+a3*la*s3*t181*t195+a2*mu*s2*t178*t181-a2*mu*s2*t178*t184-a2*mu*s2*t181*t184-a3*mu*s3*t178*t181+a3*mu*s3*t181*t184-a*a2*la*s2*s3*t195+a*a3*la*s2*s3*t195*2.0-a*a3*mu*s2*s3*t178*2.0+a*a3*mu*s2*s3*t184*2.0)/(t243+t244+t245+t246+t249+t250-a*t178*t184-a*t181*t184-s2*t178*t184*2.0-s3*t178*t184-s3*t181*t184-a*s2*s3*t184*2.0)-t258*t262*t55*t57*(a*la*t47*t51*(t259-a*s3+a2*t181-a3*t184*2.0+a*a2*s2*3.0-a*a3*s3*3.0+a2*s2*s3*2.0-a3*s2*s3)+a*t203*t260*t46*t48*(a2*s2-a3*s3))+(2.0*k*(s2*a2 - s3*a3)*a)/(b3);
    dP_dF.x2112 = t205*t46*t55*t57*(t186+t206+t207+t208+t209-t210-t211-t212-t213-t214-t215-t216+t217-t218-t219-t220-t221-t222-t223-t224-la*s2*t179*2.0-mu*s2*t179*2.0-a1*la*s2*t179*6.0+a2*la*s1*t179*6.0+a*la*t178*t180*3.0-a*la*t178*t181*3.0-a1*mu*s2*t179*3.0+a2*mu*s1*t179*3.0+a*mu*t178*t180*4.0-a*mu*t178*t181*4.0+la*s1*t178*t180-la*s1*t178*t181+la*s2*t178*t180-la*s2*t178*t181+la*s3*t178*t180*3.0-la*s3*t178*t181*3.0-la*s1*t179*t50*2.0+mu*s1*t178*t180*2.0-mu*s1*t178*t181*2.0+mu*s2*t178*t180*2.0-mu*s2*t178*t181*2.0+mu*s3*t178*t180*4.0-mu*s3*t178*t181*4.0-a*a1*la*t178*t181*7.0+a*a2*la*t178*t180*7.0+a*a2*la*t180*t181-a*a1*mu*t178*t181*4.0-a*a1*mu*t179*t181*2.0+a*a2*mu*t178*t180*4.0+a*a2*mu*t179*t180*2.0+a*a2*mu*t180*t181*2.0+a*la*s1*s3*t180-a*la*s1*s3*t181-a*la*s2*s3*t178*2.0+a*la*s2*s3*t180-a*la*s2*s3*t181+a*mu*s1*s2*t180-a*mu*s1*s2*t181+a*mu*s1*s3*t180*2.0-a*mu*s1*s3*t181*2.0+a*mu*s2*s3*t180*2.0-a*mu*s2*s3*t181*2.0-a1*la*s1*t178*t181*6.0-a1*la*s2*t178*t180*2.0+a2*la*s1*t178*t180*2.0-a1*la*s2*t178*t181*2.0+a2*la*s1*t178*t181*2.0+a2*la*s2*t178*t180*6.0-a1*la*s3*t178*t181*4.0+a2*la*s3*t178*t180*4.0+a1*la*s2*t179*t50*3.0-a2*la*s1*t179*t50*3.0-a*la*t178*t180*t50*4.0+a*la*t178*t181*t50*4.0-a1*mu*s1*t178*t181*4.0-a1*mu*s1*t179*t181*4.0-a1*mu*s2*t178*t179+a2*mu*s1*t178*t179-a1*mu*s2*t178*t180*2.0+a2*mu*s1*t178*t180*2.0-a1*mu*s2*t178*t181*2.0-a1*mu*s2*t179*t180+a2*mu*s1*t178*t181*2.0+a2*mu*s1*t179*t180-a1*mu*s2*t179*t181+a2*mu*s1*t179*t181+a1*mu*s2*t179*t184-a2*mu*s1*t179*t184-a1*mu*s2*t180*t181+a2*mu*s1*t180*t181+a2*mu*s2*t178*t180*4.0-a1*mu*s3*t178*t181*2.0+a2*mu*s2*t179*t180*4.0+a2*mu*s3*t178*t180*2.0+mu*s1*s2*s3*t180-mu*s1*s2*s3*t181-la*s1*t178*t180*t50*2.0+la*s1*t178*t181*t50*2.0-la*s2*t178*t180*t50*2.0+la*s2*t178*t181*t50*2.0-la*s3*t178*t180*t50*4.0+la*s3*t178*t181*t50*4.0-la*s1*s2*s3*t180*t50+la*s1*s2*s3*t181*t50-a*a1*la*s1*s2*t181+a*a2*la*s1*s2*t178*7.0+a*a2*la*s1*s2*t180-a*a1*la*s1*s3*t181*2.0-a*a1*la*s2*s3*t178*4.0+a*a2*la*s1*s3*t178*4.0-a*a1*la*s2*s3*t180+a*a2*la*s1*s3*t180-a*a1*la*s2*s3*t181+a*a2*la*s1*s3*t181+a*a2*la*s2*s3*t180*2.0-a*a1*mu*s1*s2*t181*2.0+a*a2*mu*s1*s2*t178*4.0+a*a2*mu*s1*s2*t179*2.0+a*a2*mu*s1*s2*t180*2.0-a*a1*mu*s2*s3*t178*2.0+a*a2*mu*s1*s3*t178*2.0-a*a1*mu*s2*s3*t180+a*a2*mu*s1*s3*t180-a*a1*mu*s2*s3*t181+a*a2*mu*s1*s3*t181+a*a1*la*t178*t181*t50*4.0+a*a1*la*t180*t181*t50*2.0-a*a2*la*t178*t180*t50*4.0+a*a1*mu*t178*t181*t184*2.0+a*a1*mu*t180*t181*t184*2.0+a*a2*mu*t178*t180*t181*2.0-a*a2*mu*t178*t180*t184*2.0+a2*la*s1*s2*s3*t178*4.0-a*la*s1*s2*t180*t50+a*la*s1*s2*t181*t50-a*la*s1*s3*t178*t50*2.0-a*la*s1*s3*t180*t50*2.0+a*la*s1*s3*t181*t50*2.0-a*la*s2*s3*t180*t50*2.0+a*la*s2*s3*t181*t50*2.0+a2*mu*s1*s2*s3*t178*2.0+a1*la*s1*t178*t181*t50*4.0+a1*la*s2*t178*t180*t50*2.0-a2*la*s1*t178*t180*t50*2.0+a1*la*s2*t178*t181*t50*2.0-a2*la*s1*t178*t181*t50*2.0+a1*la*s2*t180*t181*t50-a2*la*s1*t180*t181*t50-a2*la*s2*t178*t180*t50*4.0+a1*la*s3*t178*t181*t50*2.0-a2*la*s3*t178*t180*t50*2.0+a1*mu*s1*t178*t181*t184*4.0-a1*mu*s2*t178*t180*t181+a2*mu*s1*t178*t180*t181+a1*mu*s2*t178*t180*t184-a2*mu*s1*t178*t180*t184+a1*mu*s2*t178*t181*t184-a2*mu*s1*t178*t181*t184+a1*mu*s2*t180*t181*t184-a2*mu*s1*t180*t181*t184-a2*mu*s2*t178*t180*t184*4.0+a*a1*la*s1*s2*t178*t50*4.0+a*a1*la*s1*s2*t181*t50*2.0-a*a2*la*s1*s2*t180*t50*2.0+a*a1*la*s2*s3*t178*t50*2.0-a*a2*la*s1*s3*t178*t50*2.0+a*a1*la*s2*s3*t180*t50-a*a2*la*s1*s3*t180*t50+a*a1*la*s2*s3*t181*t50-a*a2*la*s1*s3*t181*t50-a*a1*mu*s1*s2*t178*t181*2.0+a*a1*mu*s1*s2*t178*t184*2.0+a*a1*mu*s1*s2*t181*t184*2.0+a*a2*mu*s1*s2*t178*t180*2.0-a*a2*mu*s1*s2*t180*t184*2.0+a1*la*s1*s2*s3*t178*t50*2.0)-(2.0*k*(s1*a2 - s2*a1)*a)/(b1);
    dP_dF.x3113 = -t232*(-t233+t234+t235+t251+t252+t254+t255-a*la*t180-a*mu*t180*2.0-mu*s1*t178+mu*s1*t179-mu*s1*t180+mu*s1*t184+a*a1*la*t184-a*a3*la*t180+a*a1*mu*t184-a*a3*mu*t180*2.0+a*mu*s1*s3+a1*la*s3*t178*2.0-a3*la*s1*t178*2.0+a*la*t180*t195*2.0-a*la*t184*t195*2.0+a1*mu*s3*t178*2.0-a3*mu*s1*t178*2.0+a1*mu*s3*t179-a3*mu*s1*t179+a1*mu*s3*t180-a3*mu*s1*t180+a*mu*t178*t180*2.0+la*s1*t178*t195+la*s1*t180*t195-la*s1*t184*t195-la*s3*t178*t195*2.0+mu*s1*t178*t180+mu*s3*t178*t180*2.0-a*a3*la*s1*s3+a*a1*mu*s1*s3*2.0-a*a3*mu*s1*s3-a*a1*la*t184*t195+a*a3*la*t180*t195*2.0-a*a3*mu*t178*t180*2.0+a*a3*mu*t180*t184*2.0-a*la*s1*s3*t195+a*mu*s1*s3*t178+a*mu*s1*s3*t180-a1*la*s3*t178*t195*2.0+a3*la*s1*t178*t195*2.0-a1*la*s3*t180*t195+a3*la*s1*t180*t195+a1*mu*s3*t178*t180-a3*mu*s1*t178*t180-a1*mu*s3*t178*t184+a3*mu*s1*t178*t184-a1*mu*s3*t180*t184+a3*mu*s1*t180*t184-a*a1*la*s1*s3*t195*2.0+a*a3*la*s1*s3*t195-a*a1*mu*s1*s3*t184*2.0)+t237*t242*t46*t57*(a*la*t45*t52*(t180-t184+t238-a*s3-a1*t184*2.0+a3*t180-a*a1*s3*3.0+a*a3*s1*3.0-a1*s1*s3+a3*s1*s3*2.0)-a*t203*t239*t48*t55*(a1*s3-a3*s1))-(2.0*k*(s1*a3 - s3*a1)*a)/(b2);
    dP_dF.x3223 = -(t251+t252-t253+t254+t255+t256-a*la*t181-a*mu*t181*2.0-mu*s2*t178+mu*s2*t179-mu*s2*t181+mu*s2*t184+a*a2*la*t184-a*a3*la*t181+a*a2*mu*t184-a*a3*mu*t181*2.0+a*mu*s2*s3+a2*la*s3*t178*2.0-a3*la*s2*t178*2.0+a*la*t181*t195*2.0-a*la*t184*t195*2.0+a2*mu*s3*t178*2.0-a3*mu*s2*t178*2.0+a2*mu*s3*t179-a3*mu*s2*t179+a2*mu*s3*t181-a3*mu*s2*t181+a*mu*t178*t181*2.0+la*s2*t178*t195+la*s2*t181*t195-la*s2*t184*t195-la*s3*t178*t195*2.0+mu*s2*t178*t181+mu*s3*t178*t181*2.0-a*a3*la*s2*s3+a*a2*mu*s2*s3*2.0-a*a3*mu*s2*s3-a*a2*la*t184*t195+a*a3*la*t181*t195*2.0-a*a3*mu*t178*t181*2.0+a*a3*mu*t181*t184*2.0-a*la*s2*s3*t195+a*mu*s2*s3*t178+a*mu*s2*s3*t181-a2*la*s3*t178*t195*2.0+a3*la*s2*t178*t195*2.0-a2*la*s3*t181*t195+a3*la*s2*t181*t195+a2*mu*s3*t178*t181-a3*mu*s2*t178*t181-a2*mu*s3*t178*t184+a3*mu*s2*t178*t184-a2*mu*s3*t181*t184+a3*mu*s2*t181*t184-a*a2*la*s2*s3*t195*2.0+a*a3*la*s2*s3*t195+a*a2*mu*s2*s3*t178*2.0-a*a2*mu*s2*s3*t184*2.0)/(t243+t244+t245+t246-t247-t248+t249+t250-a*t181*t184-s2*t178*t184*2.0-s3*t181*t184-a*s2*s3*t184*2.0)+t258*t262*t55*t57*(a*la*t47*t51*(t181-t184+t259-a*s3-a2*t184*2.0+a3*t181-a*a2*s3*3.0+a*a3*s2*3.0-a2*s2*s3+a3*s2*s3*2.0)-a*t203*t260*t46*t48*(a2*s3-a3*s2))-(2.0*k*(s2*a3 - s3*a2)*a)/(b3);

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
