//#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Russell Howes, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ST_VENANT_KIRCHHOFF
//##################################################################### 
#ifndef __ST_VENANT_KIRCHHOFF__
#define __ST_VENANT_KIRCHHOFF__

#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T,int d>
class ST_VENANT_KIRCHHOFF:public ISOTROPIC_CONSTITUTIVE_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;
public:
    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,d> BASE;
    using BASE::enforce_definiteness;using BASE::constant_lambda;using BASE::constant_mu;
    using BASE::constant_alpha;using BASE::constant_beta;
    using BASE::alpha;using BASE::beta;using BASE::lambda;using BASE::mu;

    T youngs_modulus,poissons_ratio;
    T failure_threshold;
    
    ST_VENANT_KIRCHHOFF(const T youngs_modulus_input=3e6,const T poissons_ratio_input=.475,const T Rayleigh_coefficient=.05);
    virtual ~ST_VENANT_KIRCHHOFF();

    T Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const override;
    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const override;
    MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const override;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int id) const override;
//#####################################################################
};
}
#endif
