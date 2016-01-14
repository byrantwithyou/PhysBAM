//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MOONEY_RIVLIN_3D2
//##################################################################### 
#ifndef __MOONEY_RIVLIN_3D2__
#define __MOONEY_RIVLIN_3D2__

#include <Tools/Math_Tools/constants.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T>
class MOONEY_RIVLIN_3D2:public ISOTROPIC_CONSTITUTIVE_MODEL<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,3> BASE;
    using BASE::enforce_definiteness;using BASE::constant_lambda;using BASE::constant_mu;
    using BASE::constant_alpha;using BASE::constant_beta;
    using BASE::Alpha;using BASE::Beta;using BASE::Lambda;using BASE::Mu;

    T mu_10,mu_01,kappa;
    T failure_threshold;
    T youngs_modulus,poissons_ratio;

    MOONEY_RIVLIN_3D2(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient);
    T Energy_Density(const DIAGONAL_MATRIX<T,3>& F,const int id) const override;
    DIAGONAL_MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const int id) const override;
    MATRIX<T,3> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& F_dot,const int id) const override;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int id) const override;

//#####################################################################
};
}
#endif
