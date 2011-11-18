//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Alexey Stomakhin, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEO_HOOKEAN_COROTATED_BLEND
//#####################################################################
#ifndef __NEO_HOOKEAN_COROTATED_BLEND__
#define __NEO_HOOKEAN_COROTATED_BLEND__

#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/CUBIC_HEAVISIDE_TRANSITION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONAL_MATRIX_DETERMINANT_HEAVISIDE_TRANSITION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/QUINTIC_HEAVISIDE_TRANSITION.h>

namespace PhysBAM{

using ::std::log;

template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d>
class NEO_HOOKEAN_COROTATED_BLEND:public ISOTROPIC_CONSTITUTIVE_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;

public:

    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,d> BASE;
    
    using BASE::enforce_definiteness;
    using BASE::constant_lambda;
    using BASE::constant_mu;
    using BASE::constant_alpha;
    using BASE::constant_beta;

    DIAGONAL_MATRIX_DETERMINANT_HEAVISIDE_TRANSITION <T,d,QUINTIC_HEAVISIDE_TRANSITION<T> > blend;

    NEO_HOOKEAN<T,d> neo_base;
    COROTATED<T,d> cor_base;

    T J_min,J_max;

public:

    NEO_HOOKEAN_COROTATED_BLEND(const T youngs_modulus=3e6,const T poissons_ratio=.475,const T Rayleigh_coefficient=.05);
    virtual ~NEO_HOOKEAN_COROTATED_BLEND();

public:

    T Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const PHYSBAM_OVERRIDE;
    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const PHYSBAM_OVERRIDE;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int triangle) const;
    void Isotropic_Stress_Derivative_Transition_Helper(DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,
        const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& neo_dP_dF,
        const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& cor_dP_dF,
        const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& DDt,
        const DIAGONAL_MATRIX<T,2>& neo_P,
        const DIAGONAL_MATRIX<T,2>& cor_P,
        const DIAGONAL_MATRIX<T,2>& Dt,
        const T neo_minus_cor,
        const T t) const;
    void Isotropic_Stress_Derivative_Transition_Helper(DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,
        const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& neo_dP_dF,
        const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& cor_dP_dF,
        const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& DDt,
        const DIAGONAL_MATRIX<T,3>& neo_P,
        const DIAGONAL_MATRIX<T,3>& cor_P,
        const DIAGONAL_MATRIX<T,3>& Dt,
        const T neo_minus_cor,
        const T t) const;
    
    MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const PHYSBAM_OVERRIDE;
    int P_From_Strain_Rate_Forces_Size() const PHYSBAM_OVERRIDE;
    void P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const PHYSBAM_OVERRIDE;
    MATRIX<T,d> P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const PHYSBAM_OVERRIDE;
};
}
#endif
