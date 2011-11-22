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
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_COROTATED_BLEND.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
NEO_HOOKEAN_COROTATED_BLEND(const T youngs_modulus,const T poissons_ratio,const T Rayleigh_coefficient):
    neo_base(youngs_modulus,poissons_ratio,Rayleigh_coefficient),
    cor_base(2*youngs_modulus,poissons_ratio,Rayleigh_coefficient),
    J_min(0.3),J_max(0.7)
{
    constant_lambda = youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    constant_mu     = youngs_modulus/(2*(1+poissons_ratio));
    constant_alpha  = Rayleigh_coefficient*constant_lambda;
    constant_beta   = Rayleigh_coefficient*constant_mu;    
    
    blend.Initialize(J_min,J_max);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
~NEO_HOOKEAN_COROTATED_BLEND()
{
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    T J = F.Determinant();

    if      (J>=J_max) return neo_base.Energy_Density(F,simplex); // Neo Hookean
    else if (J<=J_min) return cor_base.Energy_Density(F,simplex); // Corotated
    else
    {
        T t = blend.H(F);
        return t*neo_base.Energy_Density(F,simplex) + (1-t)*cor_base.Energy_Density(F,simplex);
    }
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    T J = F.Determinant();

    if      (J>=J_max) return neo_base.P_From_Strain(F,scale,simplex); // Neo Hookean
    else if (J<=J_min) return cor_base.P_From_Strain(F,scale,simplex); // Corotated
    else
    {
        T t = blend.H(F);
        return scale*(
            neo_base.P_From_Strain(F,1,simplex)*t + cor_base.P_From_Strain(F,1,simplex)*(1-t) +
            blend.DH(F)*(neo_base.Energy_Density(F,simplex) - cor_base.Energy_Density(F,simplex))
        );
    }
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int triangle) const
{
    T J = F.Determinant();

    if      (J>=J_max) neo_base.Isotropic_Stress_Derivative(F,dP_dF,triangle); // Neo Hookean
    else if (J<=J_min) cor_base.Isotropic_Stress_Derivative(F,dP_dF,triangle); // Corotated
    else
    {
        T t   = blend.H(F);
        DIAGONAL_MATRIX<T,d> Dt  = blend.DH(F);
        DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d> DDt;
        blend.DDH(F,DDt);

        T neo_minus_cor = neo_base.Energy_Density(F,triangle) - cor_base.Energy_Density(F,triangle);

        DIAGONAL_MATRIX<T,d> neo_P = neo_base.P_From_Strain(F,1,triangle);
        DIAGONAL_MATRIX<T,d> cor_P = cor_base.P_From_Strain(F,1,triangle);

        DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d> neo_dP_dF;
        DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d> cor_dP_dF;
        neo_base.Isotropic_Stress_Derivative(F,neo_dP_dF,triangle);
        cor_base.Isotropic_Stress_Derivative(F,cor_dP_dF,triangle);
        
        Isotropic_Stress_Derivative_Transition_Helper(dP_dF,neo_dP_dF,cor_dP_dF,DDt,neo_P,cor_P,Dt,neo_minus_cor,t);
    }

    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Isotropic_Stress_Derivative_Transition_Helper(DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,
    const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& neo_dP_dF,
    const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& cor_dP_dF,
    const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& DDt,
    const DIAGONAL_MATRIX<T,2>& neo_P,
    const DIAGONAL_MATRIX<T,2>& cor_P,
    const DIAGONAL_MATRIX<T,2>& Dt,
    const T neo_minus_cor,
    const T t) const
{
    dP_dF.x1111 = neo_dP_dF.x1111*t + cor_dP_dF.x1111*(1-t) + neo_minus_cor*DDt.x1111 + 2*(neo_P.x11-cor_P.x11)*Dt.x11;
    dP_dF.x2222 = neo_dP_dF.x2222*t + cor_dP_dF.x2222*(1-t) + neo_minus_cor*DDt.x2222 + 2*(neo_P.x22-cor_P.x22)*Dt.x22;
    dP_dF.x2211 = neo_dP_dF.x2211*t + cor_dP_dF.x2211*(1-t) + neo_minus_cor*DDt.x2211 + (neo_P.x22-cor_P.x22)*Dt.x11 + (neo_P.x11-cor_P.x11)*Dt.x22;
    dP_dF.x2112 = neo_dP_dF.x2112*t + cor_dP_dF.x2112*(1-t) + neo_minus_cor*DDt.x2112;
    dP_dF.x2121 = neo_dP_dF.x2121*t + cor_dP_dF.x2121*(1-t) + neo_minus_cor*DDt.x2121;
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_COROTATED_BLEND<T,d>::
Isotropic_Stress_Derivative_Transition_Helper(DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,
    const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& neo_dP_dF,
    const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& cor_dP_dF,
    const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& DDt,
    const DIAGONAL_MATRIX<T,3>& neo_P,
    const DIAGONAL_MATRIX<T,3>& cor_P,
    const DIAGONAL_MATRIX<T,3>& Dt,
    const T neo_minus_cor,
    const T t) const
{
    dP_dF.x1111 = neo_dP_dF.x1111*t + cor_dP_dF.x1111*(1-t) + neo_minus_cor*DDt.x1111 + 2*(neo_P.x11-cor_P.x11)*Dt.x11;
    dP_dF.x2222 = neo_dP_dF.x2222*t + cor_dP_dF.x2222*(1-t) + neo_minus_cor*DDt.x2222 + 2*(neo_P.x22-cor_P.x22)*Dt.x22;
    dP_dF.x3333 = neo_dP_dF.x3333*t + cor_dP_dF.x3333*(1-t) + neo_minus_cor*DDt.x3333 + 2*(neo_P.x33-cor_P.x33)*Dt.x33;

    dP_dF.x2211 = neo_dP_dF.x2211*t + cor_dP_dF.x2211*(1-t) + neo_minus_cor*DDt.x2211 + (neo_P.x22-cor_P.x22)*Dt.x11 + (neo_P.x11-cor_P.x11)*Dt.x22;
    dP_dF.x3322 = neo_dP_dF.x3322*t + cor_dP_dF.x3322*(1-t) + neo_minus_cor*DDt.x3322 + (neo_P.x33-cor_P.x33)*Dt.x22 + (neo_P.x22-cor_P.x22)*Dt.x33;
    dP_dF.x3311 = neo_dP_dF.x3311*t + cor_dP_dF.x3311*(1-t) + neo_minus_cor*DDt.x3311 + (neo_P.x33-cor_P.x33)*Dt.x11 + (neo_P.x11-cor_P.x11)*Dt.x33;

    dP_dF.x2112 = neo_dP_dF.x2112*t + cor_dP_dF.x2112*(1-t) + neo_minus_cor*DDt.x2112;
    dP_dF.x3113 = neo_dP_dF.x3113*t + cor_dP_dF.x3113*(1-t) + neo_minus_cor*DDt.x3113;
    dP_dF.x3223 = neo_dP_dF.x3223*t + cor_dP_dF.x3223*(1-t) + neo_minus_cor*DDt.x3223;

    dP_dF.x2121 = neo_dP_dF.x2121*t + cor_dP_dF.x2121*(1-t) + neo_minus_cor*DDt.x2121;
    dP_dF.x3131 = neo_dP_dF.x3131*t + cor_dP_dF.x3131*(1-t) + neo_minus_cor*DDt.x3131;
    dP_dF.x3232 = neo_dP_dF.x3232*t + cor_dP_dF.x3232*(1-t) + neo_minus_cor*DDt.x3232;
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void NEO_HOOKEAN_COROTATED_BLEND<T,d>::
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
template<class T,int d> MATRIX<T,d> NEO_HOOKEAN_COROTATED_BLEND<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
template class NEO_HOOKEAN_COROTATED_BLEND<float,2>;
template class NEO_HOOKEAN_COROTATED_BLEND<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class NEO_HOOKEAN_COROTATED_BLEND<double,2>;
template class NEO_HOOKEAN_COROTATED_BLEND<double,3>;
#endif
