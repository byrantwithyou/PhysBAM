//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Russell Howes, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED_QUARTIC.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> COROTATED_QUARTIC<T,d>::
COROTATED_QUARTIC(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T failure_threshold_input)
    :panic_threshold((T)1e-6)
{
    Set_Parameters(youngs_modulus_input,poissons_ratio_input,Rayleigh_coefficient);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> COROTATED_QUARTIC<T,d>::
~COROTATED_QUARTIC()
{
}
//#####################################################################
// Function Set_Parameters
//#####################################################################
template<class T,int d> void COROTATED_QUARTIC<T,d>::
Set_Parameters(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient)
{
    assert(poissons_ratio_input>-1&&poissons_ratio_input<.5);
    youngs_modulus=youngs_modulus_input;
    poissons_ratio=poissons_ratio_input;
    constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    constant_mu=youngs_modulus/(2*(1+poissons_ratio));
    constant_alpha=Rayleigh_coefficient*constant_lambda;
    constant_beta=Rayleigh_coefficient*constant_mu;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
// clamp to hyperbola to avoid indefiniteness "automatically"
template<class T,int d> DIAGONAL_MATRIX<T,d> COROTATED_QUARTIC<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    T c = 100;
    T scale_mu=scale*constant_mu,scale_lambda=scale*constant_lambda;
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    T Fm1Trace=Fm1.Trace();
    return 2*scale_mu*(Fm1+Fm1*Fm1*Fm1*(T)2*c)+scale_lambda*(Fm1Trace+Fm1Trace*Fm1Trace*Fm1Trace*(T)2*c);
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> COROTATED_QUARTIC<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int COROTATED_QUARTIC<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void COROTATED_QUARTIC<T,d>::
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
template<class T,int d> MATRIX<T,d> COROTATED_QUARTIC<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void COROTATED_QUARTIC<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const
{
    
    T c = 100; //this is input in by hand for now
    
    T mu=constant_mu,la=constant_lambda;//,mu2la=2*mu+la,la2mu2=2*la+2*mu;
    T t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38;
    T s1=F.x11;
    T s2=F.x22;
/*    t5 = s1+s2-2;
    t6 = s1-1;
    t7 = t5*t5;
    t8 = t7*12;
    t9 = t8+2;
    t10 = la*t9/2;
    t11 = s2-1;
    t12 = s1+s2;if(fabs(t12)<panic_threshold) t12=t12<0?-panic_threshold:panic_threshold;
    t13 = 1/t12;
    t14 = s1*s1;
    t15 = s2*s2;
    t16 = la*25;
    t17 = la*t14*2;
    t18 = la*t15*2;
    t19 = la*s1*s2*4;
    dP_dF.x1111 = t10+mu*((t6*t6)*12+2);
    dP_dF.x2222 = t10+mu*((t11*t11)*12+2);
    dP_dF.x2211 = t10;
    dP_dF.x2121 = mu*14+t16+t17+t18+t19-la*s1*12-la*s2*12-la*t13*18-mu*s1*12-mu*t13*6+mu*t14*4+mu*t15*4-mu*t13*t15*12;
    dP_dF.x2112 = -t16-t17-t18-t19+la*s1*12+la*s2*12+la*t13*18-mu*s1*12+mu*t13*6+mu*s1*s2*4+mu*t13*t14*12;*/

    t21 = s1-1;
    t22 = s1+s2-2;
    t23 = s2-1;
    t24 = t22*t22;
    t25 = c*t24*12;
    t26 = t25+2;
    t27 = la*t26/2;
    t28 = s1+s2;
    t29 = 1/t28;
    t30 = s1*s1;
    t31 = s2*s2;
    t32 = c*la*24;
    t33 = mu*2;
    t34 = c*la*12;
    t35 = c*mu*s1*4;
    t36 = t34+t35-c*la*s1*4;
    t37 = c*la*t31*2;
    t38 = c*mu*t30*4;
    dP_dF.x1111 = t27+mu*(c*(t21*t21)*12+2);
    dP_dF.x2222 = t27+mu*(c*(t23*t23)*12+2);
    dP_dF.x2211 = t27;
    dP_dF.x2121 = la+t32+t33+t37+t38+c*mu*12-la*t29*2-mu*t29*2-c*la*s1*12-c*la*s2*12-c*la*t29*16+c*la*t30*2-c*mu*s1*12-c*mu*t29*4+c*mu*t31*4+c*la*s1*s2*4-c*mu*t29*t31*12;
    dP_dF.x2112 = -la-t32-t37+t38-s1*t36+s2*t36+t29*(la*2+t33+c*la*16+c*mu*4+c*mu*t30*12)+c*la*s1*24-c*la*t30*6-c*mu*s1*12;


    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void COROTATED_QUARTIC<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const
{
    T mu=constant_mu,la=constant_lambda,mu2la=2*mu+la,la3mu2=3*la+2*mu;
    T d12=F.x11+F.x22;if(fabs(d12)<panic_threshold) d12=d12<0?-panic_threshold:panic_threshold;
    T d13=F.x11+F.x33;if(fabs(d13)<panic_threshold) d13=d13<0?-panic_threshold:panic_threshold;
    T d23=F.x22+F.x33;if(fabs(d23)<panic_threshold) d23=d23<0?-panic_threshold:panic_threshold;
    T i12=(la3mu2-la*F.x33)/d12,i13=(la3mu2-la*F.x22)/d13,i23=(la3mu2-la*F.x11)/d23;
    dPi_dF.x1111=mu2la;
    dPi_dF.x2222=mu2la;
    dPi_dF.x3333=mu2la;
    dPi_dF.x2211=la;
    dPi_dF.x3311=la;
    dPi_dF.x3322=la;
    dPi_dF.x2112=i12-la;
    dPi_dF.x3113=i13-la;
    dPi_dF.x3223=i23-la;
    dPi_dF.x2121=mu2la-i12;
    dPi_dF.x3131=mu2la-i13;
    dPi_dF.x3232=mu2la-i23;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T COROTATED_QUARTIC<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    T c = 100;
    DIAGONAL_MATRIX<T,d> Fm1=F-1;
    return constant_mu*(Fm1*Fm1+Fm1*Fm1*Fm1*Fm1*c).Trace()+(T).5*constant_lambda*(sqr(Fm1.Trace())*sqr(Fm1.Trace())*c+sqr(Fm1.Trace()));
}
template class COROTATED_QUARTIC<float,2>;
template class COROTATED_QUARTIC<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COROTATED_QUARTIC<double,2>;
template class COROTATED_QUARTIC<double,3>;
#endif
