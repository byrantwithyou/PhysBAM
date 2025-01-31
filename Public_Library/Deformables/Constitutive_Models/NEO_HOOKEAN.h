//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEO_HOOKEAN
//#####################################################################
#ifndef __NEO_HOOKEAN__
#define __NEO_HOOKEAN__

#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

using ::std::log;

template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d>
class NEO_HOOKEAN:public ISOTROPIC_CONSTITUTIVE_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;
public:
    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,d> BASE;
    using BASE::enforce_definiteness;using BASE::constant_lambda;using BASE::constant_mu;
    using BASE::constant_alpha;using BASE::constant_beta;
    using BASE::Alpha;using BASE::Beta;using BASE::Lambda;using BASE::Mu;

    T youngs_modulus,poissons_ratio;
    T failure_threshold;
    bool use_constant_ife;
private:
    T dth_root_failure_threshold;
public:

    NEO_HOOKEAN(const T youngs_modulus_input=3e6,const T poissons_ratio_input=.475,const T Rayleigh_coefficient=.05,const T failure_threshold_input=.25);
    virtual ~NEO_HOOKEAN();

private:
    DIAGONAL_MATRIX<T,1> Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,1>& F) const;
    DIAGONAL_MATRIX<T,2> Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,2>& F) const;
    DIAGONAL_MATRIX<T,3> Clamp_To_Hyperbola(const DIAGONAL_MATRIX<T,3>& F) const;
public:

    // clamp to hyperbola to avoid indefiniteness "automatically"
    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const override;
    T Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const override;

    MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const override;
    void Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,1>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<T,1> >& dP_dF,const int id) const;
    void Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<T,2> >& dP_dF,const int id) const;
    void Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<T,3> >& dPi_dF,const int id) const;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dPi_dF,const int id) const override;
    int P_From_Strain_Rate_Forces_Size() const override;
    void P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const int id) const override;
    MATRIX<T,d> P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const int id) const override;
    virtual T Robust_Divided_Pressure(T J,const int id) const override;
    virtual T Pressure_Bound(T J,const int id) const override;
//#####################################################################
};
}
#endif
