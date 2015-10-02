//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COROTATED_FIXED
//#####################################################################
#ifndef __COROTATED_FIXED__
#define __COROTATED_FIXED__

#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

using ::std::log;

template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d>
class COROTATED_FIXED:public ISOTROPIC_CONSTITUTIVE_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;
public:
    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,d> BASE;
    using BASE::enforce_definiteness;using BASE::constant_lambda;using BASE::constant_mu;using BASE::constant_alpha;using BASE::constant_beta;

    T youngs_modulus,poissons_ratio;
    T panic_threshold;

    COROTATED_FIXED(const T youngs_modulus_input=3e6,const T poissons_ratio_input=.475,const T Rayleigh_coefficient=.05,const T failure_threshold_input=.25);
    virtual ~COROTATED_FIXED();

    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const override;
    T Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const override;
    void Set_Parameters(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient);
    void Zero_Out_Mu();

    MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const override;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,1>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,1>& dP_dF,const int triangle) const;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int triangle) const;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron) const;
    int P_From_Strain_Rate_Forces_Size() const override;
    void P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const override;
    MATRIX<T,d> P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const override;

//#####################################################################
};
}
#endif
