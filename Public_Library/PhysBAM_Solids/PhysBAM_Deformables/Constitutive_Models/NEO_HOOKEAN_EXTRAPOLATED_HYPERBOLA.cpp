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
NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input):
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
template<class T,int d> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
~NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA()
{
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
    
//    T t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t0;

    T t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t0;
    
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
    t17 = 1.0/(t9*t9);
    t18 = t16*t16;
    t19 = 1.0/(t12*t12);
    t0 = -mu*t14+t16*(mu*t12-mu*t15+la*t14*t15)+t16*(-mu*t11+mu*t9+la*t11*t14)+la*(t14*t14)*(1.0/2.0)+mu*(t12*t12+t9*t9-2.0)*(1.0/2.0)+t18*(mu+la*t17+mu*t17-la*t14*t17)*(1.0/2.0)+t18*(mu+la*t19+mu*t19-la*t14*t19)*(1.0/2.0)+la*t11*t15*t18;

/*    t40 = s1-s2;
    t41 = s2*(1.0/2.0);
    t42 = c*4.0;
    t43 = t40*t40;
    t44 = t42+t43;
    t45 = sqrt(t44);
    t46 = t45*(1.0/2.0);
    t48 = s1*(1.0/2.0);
    t47 = t41+t46-t48;
    t49 = 1.0/t47;
    t50 = -t41+t46+t48;
    t51 = t47*t50;
    t52 = log(t51);
    t53 = 1.0/t50;
    t54 = t41-t46+t48;
    t55 = 1.0/(t47*t47);
    t56 = t54*t54;
    t57 = 1.0/(t50*t50);
    t0 = -mu*t52-t54*(mu*t47-mu*t49+la*t49*t52)-t54*(mu*t50-mu*t53+la*t52*t53)+la*(t52*t52)*(1.0/2.0)+mu*(t47*t47+t50*t50-2.0)*(1.0/2.0)+t56*(mu+la*t55+mu*t55-la*t52*t55)*(1.0/2.0)+t56*(mu+la*t57+mu*t57-la*t52*t57)*(1.0/2.0)+la*t49*t53*t56;*/
    
    return t0;

}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return 0;
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
    T t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t80,t81,t82,t83,t84,t85,t86,t87;
   /* T t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t118,t119,t120,t121,t122,t123,t124,t125,t126,t127,t128,t129,t130,t131,t132,t133,t134,t135,t136,t137,t138,t139,t140,t141,t142,t143,t144,t145,t146,t147,t148,t149,t150,t151,t152,t153;
    
    
    

        t95 = s1-s2;
        t96 = c*4.0;
        t97 = t95*t95;
        t98 = t96+t97;
        t99 = s1*2.0;
        t100 = s2*2.0;
        t101 = 1.0/sqrt(t98);
        t102 = t100-t99;
        t103 = t101*t102*(1.0/4.0);
        t104 = s1*(1.0/2.0);
        t105 = s2*(1.0/2.0);
        t106 = sqrt(t98);
        t107 = t103+1.0/2.0;
        t108 = t106*(1.0/2.0);
        t109 = -t104+t105+t108;
        t110 = 1.0/(t109*t109);
        t111 = t104-t105+t108;
        t112 = t103-1.0/2.0;
        t113 = t107*t111;
        t114 = t109*t112;
        t115 = 1.0/(t111*t111);
        t116 = t113+t114;
        t117 = t109*t111;
        t118 = log(t117);
        t119 = t101*t102*(1.0/2.0);
        t120 = 1.0/t109;
        t121 = 1.0/t111;
        t122 = t104+t105-t108;
        t123 = 1.0/(t109*t109*t109);
        t124 = t122*t122;
        t125 = 1.0/(t111*t111*t111);
        t126 = mu*t107;
        t127 = mu*t107*t110;
        t128 = la*t110*t116*t121;
        t129 = mu*t112;
        t130 = mu*t112*t115;
        t131 = la*t115*t116*t120;
        t132 = t119+1.0;
        t133 = t109*t132;
        t134 = t119-1.0;
        t135 = t111*t134;
        t136 = mu*t109;
        t137 = la*t118*t120;
        t138 = t136+t137-mu*t120;
        t139 = mu*t111;
        t140 = la*t118*t121;
        t141 = t139+t140-mu*t121;
        t142 = la*t107*t123*2.0;
        t143 = mu*t107*t123*2.0;
        t144 = la*t116*t121*t123;
        t145 = la*t112*t125*2.0;
        t146 = mu*t112*t125*2.0;
        t147 = la*t116*t120*t125;
        t148 = la*t110;
        t149 = mu*t110;
        t150 = mu+t148+t149-la*t110*t118;
        t151 = la*t115;
        t152 = mu*t115;
        t153 = mu+t151+t152-la*t115*t118;
        result.x11 = mu*(t133+t135)*(-1.0/2.0)+t122*(t129+t130+t131-la*t115*t118*(t103-1.0/2.0))+t124*(t142+t143+t144-la*t118*t123*(t103+1.0/2.0)*2.0)*(1.0/2.0)+t124*(t145+t146+t147-la*t118*t125*(t103-1.0/2.0)*2.0)*(1.0/2.0)-t107*t138-t107*t141+(t104+t105-t106*(1.0/2.0))*(t126+t127+t128-la*t110*t118*(t103+1.0/2.0))+t122*t150*(t103+1.0/2.0)+t122*t153*(t103+1.0/2.0)+mu*t120*t121*(t113+t114)-la*t116*t118*t120*t121+la*t110*t121*t124*(t103+1.0/2.0)+la*t115*t120*t124*(t103-1.0/2.0)+la*t120*t121*t122*(t103+1.0/2.0)*2.0;
        result.x22 = mu*(t133+t135)*(1.0/2.0)+t138*(t103-1.0/2.0)+t141*(t103-1.0/2.0)-t122*(t126+t127+t128-la*t107*t110*t118)-t122*(t129+t130+t131-la*t112*t115*t118)-t124*(t142+t143+t144-la*t107*t118*t123*2.0)*(1.0/2.0)-t124*(t145+t146+t147-la*t112*t118*t125*2.0)*(1.0/2.0)-t112*t122*t150-t112*t122*t153-mu*t116*t120*t121+la*t118*t120*t121*(t113+t114)-la*t107*t110*t121*t124-la*t112*t115*t120*t124-la*t112*t120*t121*t122*2.0;*/
    
    t21 = s1-s2;
    t22 = c*4.0;
    t23 = t21*t21;
    t24 = t22+t23;
    t25 = s1*2.0;
    t26 = s2*2.0;
    t27 = t25-t26;
    t28 = 1.0/sqrt(t24);
    t29 = t27*t28*(1.0/4.0);
    t30 = t29-1.0/2.0;
    t31 = s1*(1.0/2.0);
    t32 = s2*(1.0/2.0);
    t33 = sqrt(t24);
    t34 = t33*(1.0/2.0);
    t35 = -t31+t32+t34;
    t36 = 1.0/(t35*t35);
    t37 = t31-t32+t34;
    t38 = t29+1.0/2.0;
    t39 = t30*t37;
    t40 = t35*t38;
    t41 = t39+t40;
    t42 = 1.0/(t37*t37);
    t43 = t35*t37;
    t44 = log(t43);
    t45 = t27*t28*(1.0/2.0);
    t46 = 1.0/t35;
    t47 = 1.0/t37;
    t48 = t31+t32-t34;
    t49 = 1.0/(t35*t35*t35);
    t50 = t48*t48;
    t51 = 1.0/(t37*t37*t37);
    t52 = mu*t30;
    t53 = mu*t30*t36;
    t54 = la*t36*t41*t47;
    t55 = t52+t53+t54-la*t30*t36*t44;
    t56 = mu*t38;
    t57 = mu*t38*t42;
    t58 = la*t41*t42*t46;
    t59 = t56+t57+t58-la*t38*t42*t44;
    t60 = t48*t59;
    t61 = t45-1.0;
    t62 = t35*t61;
    t63 = t45+1.0;
    t64 = t37*t63;
    t65 = t62+t64;
    t66 = mu*t65*(1.0/2.0);
    t67 = mu*t35;
    t68 = la*t44*t46;
    t69 = t67+t68-mu*t46;
    t70 = mu*t37;
    t71 = la*t44*t47;
    t72 = t70+t71-mu*t47;
    t73 = la*t30*t49*2.0;
    t74 = mu*t30*t49*2.0;
    t75 = la*t41*t47*t49;
    t76 = t73+t74+t75-la*t30*t44*t49*2.0;
    t77 = la*t38*t51*2.0;
    t78 = mu*t38*t51*2.0;
    t79 = la*t41*t46*t51;
    t80 = t77+t78+t79-la*t38*t44*t51*2.0;
    t81 = la*t36;
    t82 = mu*t36;
    t83 = mu+t81+t82-la*t36*t44;
    t84 = la*t42;
    t85 = mu*t42;
    t86 = mu+t84+t85-la*t42*t44;
    t87 = la*t41*t44*t46*t47;
    result.x11 = t60+t66+t87-t30*t69-t30*t72-t50*t76*(1.0/2.0)-t50*t80*(1.0/2.0)+t55*(t31+t32-t33*(1.0/2.0))-t30*t48*t83-t30*t48*t86-mu*t41*t46*t47-la*t30*t36*t47*t50-la*t30*t46*t47*t48*2.0-la*t38*t42*t46*t50;
    result.x22 = -t60-t66-t87-t48*t55+t38*t69+t38*t72+t50*t76*(1.0/2.0)+t50*t80*(1.0/2.0)+t38*t48*t83+t38*t48*t86+mu*t41*t46*t47+la*t30*t36*t47*t50+la*t38*t46*t47*t48*2.0+la*t38*t42*t46*t50;


    return scale*result;

}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const
{
    PHYSBAM_FATAL_ERROR();
    return DIAGONAL_MATRIX<T,3>();
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
   /* T t155,t156,t157,t158,t159,t160,t161,t162,t163,t164,t165,t166,t167,t168,t169,t170,t171,t172,t173,t174,t175,t176,t177,t178,t179,t180,t181,t182,t183,t184,t185,t186,t187;
    
    t155 = 1.0/c;
    t156 = s1*s1;
    t157 = s2*s2;
    t158 = c*4.0;
    t162 = s1*s2*2.0;
    t159 = t156+t157+t158-t162;
    t160 = 1.0/pow(t159,3.0/2.0);
    t161 = 1.0/(c*c);
    t163 = log(c);
    t164 = t156*t156;
    t165 = t157*t157;
    t166 = mu*3.0;
    t167 = la*t155*3.0;
    t168 = mu*t155;
    t169 = mu*s2*t157*t160*2.0;
    t170 = mu*s2*t156*t160*6.0;
    t171 = la*s1*s2*t161*t163*3.0;
    t172 = la*s1*t160*1.2E1;
    t173 = la*s2*t160*1.2E1;
    t174 = mu*s1*t160*8.0;
    t175 = mu*s2*t160*8.0;
    t176 = c*mu*s1*t160*4.0;
    t177 = c*mu*s2*t160*4.0;
    t178 = c*c;
    t179 = sqrt(t159);
    t180 = la*s1*t156;
    t181 = la*s2*t157;
    t182 = mu*s1*t156;
    t183 = mu*s2*t157;
    t184 = la*t164;
    t185 = la*t165;
    t186 = mu*t164;
    t187 = mu*t165;
    dP_dF.x1111 = t166+t167+t168+t169+t170+t171+t173+t175+t177-la*s1*t160*3.6E1+la*t156*t161*3.0-la*t155*t163+la*t157*t161-mu*s1*t160*2.4E1+mu*t156*t161*3.0+mu*t157*t161-c*mu*s1*t160*1.2E1-la*s1*s2*t161*3.0+la*s1*t160*t163*2.4E1-la*s2*t160*t163*8.0-mu*s1*s2*t161*3.0-la*t156*t161*t163*3.0-la*t157*t161*t163-mu*s1*t156*t160*2.0-mu*s1*t157*t160*6.0-la*s1*t155*t156*t160*2.1E1-la*s1*t155*t157*t160*2.7E1-la*s1*t160*t161*t164*3.0-la*s1*t160*t161*t165*3.0+la*s2*t155*t156*t160*4.5E1+la*s2*t155*t157*t160*3.0+la*s2*t160*t161*t164*1.2E1-mu*s1*t155*t156*t160*1.9E1-mu*s1*t155*t157*t160*2.1E1-mu*s1*t160*t161*t164*3.0-mu*s1*t160*t161*t165*3.0+mu*s2*t155*t156*t160*3.9E1+mu*s2*t155*t157*t160+mu*s2*t160*t161*t164*1.2E1+la*s1*t155*t156*t160*t163*1.9E1-la*s1*t156*t157*t160*t161*1.8E1+la*s1*t155*t157*t160*t163*2.1E1+la*s1*t160*t161*t163*t164*3.0+la*s1*t160*t161*t163*t165*3.0-la*s2*t155*t156*t160*t163*3.9E1+la*s2*t156*t157*t160*t161*1.2E1-la*s2*t155*t157*t160*t163-la*s2*t160*t161*t163*t164*1.2E1-mu*s1*t156*t157*t160*t161*1.8E1+mu*s2*t156*t157*t160*t161*1.2E1+la*s1*t156*t157*t160*t161*t163*1.8E1-la*s2*t156*t157*t160*t161*t163*1.2E1;
    dP_dF.x2222 = t166+t167+t168-t169-t170+t171+t172+t174+t176-la*s2*t160*3.6E1+la*t156*t161-la*t155*t163+la*t157*t161*3.0-mu*s2*t160*2.4E1+mu*t156*t161+mu*t157*t161*3.0-c*mu*s2*t160*1.2E1-la*s1*s2*t161*3.0-la*s1*t160*t163*8.0+la*s2*t160*t163*2.4E1-mu*s1*s2*t161*3.0-la*t156*t161*t163-la*t157*t161*t163*3.0+mu*s1*t156*t160*2.0+mu*s1*t157*t160*6.0+la*s1*t155*t156*t160*3.0+la*s1*t155*t157*t160*4.5E1+la*s1*t160*t161*t165*1.2E1-la*s2*t155*t156*t160*2.7E1-la*s2*t155*t157*t160*2.1E1-la*s2*t160*t161*t164*3.0-la*s2*t160*t161*t165*3.0+mu*s1*t155*t156*t160+mu*s1*t155*t157*t160*3.9E1+mu*s1*t160*t161*t165*1.2E1-mu*s2*t155*t156*t160*2.1E1-mu*s2*t155*t157*t160*1.9E1-mu*s2*t160*t161*t164*3.0-mu*s2*t160*t161*t165*3.0-la*s1*t155*t156*t160*t163+la*s1*t156*t157*t160*t161*1.2E1-la*s1*t155*t157*t160*t163*3.9E1-la*s1*t160*t161*t163*t165*1.2E1+la*s2*t155*t156*t160*t163*2.1E1-la*s2*t156*t157*t160*t161*1.8E1+la*s2*t155*t157*t160*t163*1.9E1+la*s2*t160*t161*t163*t164*3.0+la*s2*t160*t161*t163*t165*3.0+mu*s1*t156*t157*t160*t161*1.2E1-mu*s2*t156*t157*t160*t161*1.8E1-la*s1*t156*t157*t160*t161*t163*1.2E1+la*s2*t156*t157*t160*t161*t163*1.8E1;
    dP_dF.x2211 = mu*-2.0+t172+t173+t174+t175+t176+t177-la*t155-la*t156*t161*(3.0/2.0)-la*t157*t161*(3.0/2.0)-mu*t156*t161*(3.0/2.0)-mu*t157*t161*(3.0/2.0)+la*s1*s2*t161*2.0-la*s1*t160*t163*8.0-la*s2*t160*t163*8.0+mu*s1*s2*t161*2.0+la*t156*t161*t163*(3.0/2.0)+la*t157*t161*t163*(3.0/2.0)-la*s1*s2*t161*t163*2.0+la*s1*t155*t156*t160*9.0-la*s1*t155*t157*t160*9.0+la*s1*t160*t161*t164*(3.0/2.0)-la*s1*t160*t161*t165*(9.0/2.0)-la*s2*t155*t156*t160*9.0+la*s2*t155*t157*t160*9.0-la*s2*t160*t161*t164*(9.0/2.0)+la*s2*t160*t161*t165*(3.0/2.0)+mu*s1*t155*t156*t160*9.0-mu*s1*t155*t157*t160*9.0+mu*s1*t160*t161*t164*(3.0/2.0)-mu*s1*t160*t161*t165*(9.0/2.0)-mu*s2*t155*t156*t160*9.0+mu*s2*t155*t157*t160*9.0-mu*s2*t160*t161*t164*(9.0/2.0)+mu*s2*t160*t161*t165*(3.0/2.0)-la*s1*t155*t156*t160*t163*9.0+la*s1*t156*t157*t160*t161*3.0+la*s1*t155*t157*t160*t163*9.0-la*s1*t160*t161*t163*t164*(3.0/2.0)+la*s1*t160*t161*t163*t165*(9.0/2.0)+la*s2*t155*t156*t160*t163*9.0+la*s2*t156*t157*t160*t161*3.0-la*s2*t155*t157*t160*t163*9.0+la*s2*t160*t161*t163*t164*(9.0/2.0)-la*s2*t160*t161*t163*t165*(3.0/2.0)+mu*s1*t156*t157*t160*t161*3.0+mu*s2*t156*t157*t160*t161*3.0-la*s1*t156*t157*t160*t161*t163*3.0-la*s2*t156*t157*t160*t161*t163*3.0;
    dP_dF.x2121 = (t180+t181+t182+t183+c*la*s1*3.0+c*la*s2*3.0+c*mu*s1+c*mu*s2+mu*s1*t178*3.0+mu*s2*t178*3.0-c*la*s1*t163-c*la*s2*t163-la*s1*t156*t163-la*s2*t157*t163)/(s1*t178+s2*t178)-(t184+t185+t186+t187+la*t178*4.0+c*la*t156*5.0+c*la*t157*5.0+c*mu*t156*3.0+c*mu*t157*3.0+c*mu*t178*4.0-la*t163*t164-la*t163*t165+mu*t156*t178*2.0+mu*t157*t178*2.0+c*la*s1*s2*2.0+c*mu*s1*s2*2.0-c*la*t156*t163*3.0-c*la*t157*t163*3.0-la*s1*s2*t156-la*s1*s2*t157-mu*s1*s2*t156-mu*s1*s2*t157-c*la*s1*s2*t163*2.0+la*s1*s2*t156*t163+la*s1*s2*t157*t163)/(s1*t178*t179+s2*t178*t179);
    dP_dF.x2112 = -(t184+t185+t186+t187-la*t178*8.0+c*la*t156*2.0+c*la*t157*2.0+c*mu*t156*2.0+c*mu*t157*2.0-c*mu*t178*8.0-la*t156*t157*6.0-la*t163*t164-la*t163*t165-mu*t156*t157*6.0+c*la*s1*s2*2.0E1+c*mu*s1*s2*1.2E1-c*la*t156*t163*2.0-c*la*t157*t163*2.0+la*s1*s2*t156*2.0+la*s1*s2*t157*2.0+mu*s1*s2*t156*2.0+mu*s1*s2*t157*2.0+mu*s1*s2*t178*8.0+la*t156*t157*t163*6.0-c*la*s1*s2*t163*1.2E1-la*s1*s2*t156*t163*2.0-la*s1*s2*t157*t163*2.0)/(s1*t178*t179*2.0+s2*t178*t179*2.0)+(t180+t181+t182+t183+c*la*s1*2.0+c*la*s2*2.0+la*s1*t157+la*s2*t156+mu*s1*t157+mu*s1*t178*4.0+mu*s2*t156+mu*s2*t178*4.0-la*s1*t156*t163-la*s1*t157*t163-la*s2*t156*t163-la*s2*t157*t163)/(s1*t178*2.0+s2*t178*2.0);*/
    
    T t89,t90,t91,t92,t93,t94,t95,t96,t97,t98,t99,t100,t101,t102,t103,t104,t105,t106,t107;
    
    t89 = 1.0/c;
    t90 = 1.0/(c*c);
    t91 = s1*s1;
    t92 = s2*s2;
    t93 = log(c);
    t94 = c*4.0;
    t97 = s1*s2*2.0;
    t95 = t91+t92+t94-t97;
    t96 = 1.0/sqrt(t95);
    t98 = la*t89*3.0;
    t99 = mu*t89*3.0;
    t100 = la*s1*s2*t90*t93*3.0;
    t101 = la*s1*t89*t96*3.0;
    t102 = la*s2*t89*t96*3.0;
    t103 = mu*s1*t89*t96*3.0;
    t104 = mu*s2*t89*t96*3.0;
    t105 = pow(t95,3.0/2.0);
    t106 = c*c;
    t107 = sqrt(t95);
    dP_dF.x1111= mu+t100+t102+t104+t98+t99-la*t89*t93*3.0+la*t90*t91*3.0+la*t90*t92+mu*t90*t91*3.0+mu*t90*t92-la*s1*s2*t90*3.0-la*s1*t89*t96*9.0-mu*s1*s2*t90*3.0-la*t90*t91*t93*3.0-la*t90*t92*t93-mu*s1*t89*t96*9.0+la*s1*t89*t93*t96*9.0-la*s1*t90*t91*t96*3.0-la*s1*t90*t92*t96*3.0-la*s2*t89*t93*t96*3.0+la*s2*t90*t91*t96*6.0-mu*s1*t90*t91*t96*3.0-mu*s1*t90*t92*t96*3.0+mu*s2*t90*t91*t96*6.0+la*s1*t90*t91*t93*t96*3.0+la*s1*t90*t92*t93*t96*3.0-la*s2*t90*t91*t93*t96*6.0;
    dP_dF.x2222 = mu+t100+t101+t103+t98+t99-la*t89*t93*3.0+la*t90*t91+la*t90*t92*3.0+mu*t90*t91+mu*t90*t92*3.0-la*s1*s2*t90*3.0-la*s2*t89*t96*9.0-mu*s1*s2*t90*3.0-la*t90*t91*t93-la*t90*t92*t93*3.0-mu*s2*t89*t96*9.0-la*s1*t89*t93*t96*3.0+la*s1*t90*t92*t96*6.0+la*s2*t89*t93*t96*9.0-la*s2*t90*t91*t96*3.0-la*s2*t90*t92*t96*3.0+mu*s1*t90*t92*t96*6.0-mu*s2*t90*t91*t96*3.0-mu*s2*t90*t92*t96*3.0-la*s1*t90*t92*t93*t96*6.0+la*s2*t90*t91*t93*t96*3.0+la*s2*t90*t92*t93*t96*3.0;
    dP_dF.x2211 = t101+t102+t103+t104-la*t89-mu*t89*2.0+la*t89*t93*2.0-la*t90*t91*(3.0/2.0)-la*t90*t92*(3.0/2.0)-mu*t90*t91*(3.0/2.0)-mu*t90*t92*(3.0/2.0)+la*s1*s2*t90*2.0+mu*s1*s2*t90*2.0+la*t90*t91*t93*(3.0/2.0)+la*t90*t92*t93*(3.0/2.0)-la*s1*s2*t90*t93*2.0-la*s1*t89*t93*t96*3.0+la*s1*t90*t91*t96*(3.0/2.0)-la*s1*t90*t92*t96*(3.0/2.0)-la*s2*t89*t93*t96*3.0-la*s2*t90*t91*t96*(3.0/2.0)+la*s2*t90*t92*t96*(3.0/2.0)+mu*s1*t90*t91*t96*(3.0/2.0)-mu*s1*t90*t92*t96*(3.0/2.0)-mu*s2*t90*t91*t96*(3.0/2.0)+mu*s2*t90*t92*t96*(3.0/2.0)-la*s1*t90*t91*t93*t96*(3.0/2.0)+la*s1*t90*t92*t93*t96*(3.0/2.0)+la*s2*t90*t91*t93*t96*(3.0/2.0)-la*s2*t90*t92*t93*t96*(3.0/2.0);
    dP_dF.x2121 = (t90*(-la*t105-mu*t105+c*la*s1*1.2E1+c*la*s2*1.2E1+c*mu*s1*1.2E1+c*mu*s2*1.2E1+la*s1*t91*4.0+la*s2*t92*4.0-la*t107*t91*3.0-la*t107*t92*3.0+la*t105*t93+mu*s1*t106*4.0+mu*s2*t106*4.0+mu*s1*t91*4.0+mu*s2*t92*4.0-mu*t107*t91*3.0-mu*t107*t92*3.0-c*la*s1*t93*1.2E1-c*la*s2*t93*1.2E1-la*s1*s2*t107*6.0-la*s1*t91*t93*4.0-la*s2*t92*t93*4.0-mu*s1*s2*t107*6.0+la*t107*t91*t93*3.0+la*t107*t92*t93*3.0+la*s1*s2*t107*t93*6.0)*(1.0/4.0))/(s1+s2);
    dP_dF.x2112 = (la*t105*(1.0/4.0)+mu*t105*(1.0/4.0)+c*(la*s1+la*s2+mu*s1*2.0+mu*s2*2.0-la*s1*t93*2.0-la*s2*t93*2.0)+la*s1*t91*(1.0/2.0)+la*s1*t92*(1.0/2.0)+la*s2*t91*(1.0/2.0)+la*s2*t92*(1.0/2.0)-la*t107*t91*(3.0/4.0)-la*t107*t92*(3.0/4.0)-la*t105*t93*(1.0/4.0)+mu*s1*t91*(1.0/2.0)+mu*s1*t92*(1.0/2.0)+mu*s2*t91*(1.0/2.0)+mu*s2*t92*(1.0/2.0)-mu*t107*t91*(3.0/4.0)-mu*t107*t92*(3.0/4.0)-la*s1*s2*t107*(3.0/2.0)-la*s1*t91*t93*(1.0/2.0)-la*s1*t92*t93*(1.0/2.0)-la*s2*t91*t93*(1.0/2.0)-la*s2*t92*t93*(1.0/2.0)-mu*s1*s2*t107*(3.0/2.0)+la*t107*t91*t93*(3.0/4.0)+la*t107*t92*t93*(3.0/4.0)+la*s1*s2*t107*t93*(3.0/2.0))/(s1*t106+s2*t106);


}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int triangle) const
{
    PHYSBAM_FATAL_ERROR();
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
