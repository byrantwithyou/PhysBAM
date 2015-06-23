//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Alexey Stomakhin, Russell Howes, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEO_HOOKEAN_EXTRAPOLATED2
//#####################################################################
#ifndef __NEO_HOOKEAN_EXTRAPOLATED2__
#define __NEO_HOOKEAN_EXTRAPOLATED2__

#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <Deformables/Constitutive_Models/NEO_HOOKEAN_ENERGY.h>
namespace PhysBAM{

using ::std::log;

template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d>
class NEO_HOOKEAN_EXTRAPOLATED2:public ISOTROPIC_CONSTITUTIVE_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;

public:

    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,d> BASE;
    
    using BASE::enforce_definiteness;
    using BASE::constant_lambda;
    using BASE::constant_mu;
    using BASE::constant_alpha;
    using BASE::constant_beta;

    NEO_HOOKEAN_ENERGY<T> base;

    T youngs_modulus,poissons_ratio;
    T extrapolation_cutoff;
    T extra_force_coefficient;
    T panic_threshold;

public:

    NEO_HOOKEAN_EXTRAPOLATED2(const T youngs_modulus_input=3e6,const T poissons_ratio_input=.475,const T Rayleigh_coefficient=.05,const T extrapolation_cutoff=.3,const T extra_force_coefficient=20);
    virtual ~NEO_HOOKEAN_EXTRAPOLATED2();

public:

    void Update_Lame_Constants(const T youngs_modulus_input, const T poissons_ratio_input,const T Rayleigh_coefficient_input) override;
    T Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const override;
    T Energy_Density_Helper(const DIAGONAL_MATRIX<T,2>& F,const int simplex) const;
    T Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int simplex) const;
    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const override;
    DIAGONAL_MATRIX<T,2> P_From_Strain_Helper(const DIAGONAL_MATRIX<T,2>& F,const T scale,const int simplex) const;
    DIAGONAL_MATRIX<T,3> P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int triangle) const;
    void Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const;
    void Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int triangle) const;
    
    MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const override;
    int P_From_Strain_Rate_Forces_Size() const override;
    void P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const override;
    MATRIX<T,d> P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const override;
};
}
#endif
