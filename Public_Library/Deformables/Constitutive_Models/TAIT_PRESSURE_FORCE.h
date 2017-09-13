//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TAIT_PRESSURE_FORCE
//#####################################################################
#ifndef __TAIT_PRESSURE_FORCE__
#define __TAIT_PRESSURE_FORCE__

#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

using ::std::log;

template<class T,int d> class DIAGONAL_MATRIX;
template<class TV>
class TAIT_PRESSURE_FORCE:public ISOTROPIC_CONSTITUTIVE_MODEL<typename TV::SCALAR,TV::m>
{
    typedef typename TV::SCALAR T;
public:
    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,TV::m> BASE;
    using BASE::enforce_definiteness;using BASE::constant_lambda;using BASE::constant_mu;
    using BASE::constant_alpha;using BASE::constant_beta;
    using BASE::Alpha;using BASE::Beta;using BASE::Lambda;using BASE::Mu;

    T stiffness,tait_const;

    TAIT_PRESSURE_FORCE(const T stiffness_input,const T tait_const_input);
    virtual ~TAIT_PRESSURE_FORCE();

    // clamp to hyperbola to avoid indefiniteness "automatically"
    DIAGONAL_MATRIX<T,TV::m> P_From_Strain(const DIAGONAL_MATRIX<T,TV::m>& F,const int id) const override;
    T Energy_Density(const DIAGONAL_MATRIX<T,TV::m>& F,const int id) const override;

    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,TV::m>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dP_dF,const int id) const override;
//#####################################################################
};
}
#endif
