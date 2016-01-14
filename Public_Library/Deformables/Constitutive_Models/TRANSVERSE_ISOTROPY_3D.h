//#####################################################################
// Copyright 2005, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRANSVERSE_ISOTROPY_3D
//##################################################################### 
#ifndef __TRANSVERSE_ISOTROPY_3D__
#define __TRANSVERSE_ISOTROPY_3D__

#include <Tools/Arrays/ARRAY.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Deformables/Constitutive_Models/ANISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T>
class TRANSVERSE_ISOTROPY_3D:public ANISOTROPIC_CONSTITUTIVE_MODEL<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    typedef ANISOTROPIC_CONSTITUTIVE_MODEL<T,3> BASE;
    using BASE::enforce_definiteness;using BASE::constant_alpha;using BASE::constant_beta;
    using BASE::use_isotropic_component_of_stress_derivative_only;
    using BASE::Alpha;using BASE::Beta;using BASE::Lambda;using BASE::Mu;
    
    T failure_threshold;
    ARRAY<VECTOR<T,3> > fiber_field;

    TRANSVERSE_ISOTROPY_3D(const T failure_threshold_input);
    void Initialize_Fiber_Field_From_Current_State(const STRAIN_MEASURE<TV,3>& strain_measure,const ARRAY<VECTOR<T,3> >& fiber_field_input);
    inline int Hessian_Index(const int m,const int n) const;
    void Invariants(ARRAY<T>& invariants,const DIAGONAL_MATRIX<T,3>& C,const VECTOR<T,3> V_fiber) const;
    SYMMETRIC_MATRIX<T,3> S(const DIAGONAL_MATRIX<T,3>& C,const VECTOR<T,3>& V_fiber,const ARRAY<T>& invariants,const ARRAY<T> energy_gradient) const;
    MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const int id) const override;
    MATRIX<T,3> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& F_dot,const int id) const override;
    void Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,DIAGONALIZED_STRESS_DERIVATIVE<T,3>& dP_dF,const int id) const override;

//#####################################################################
    virtual void Energy_Gradient(ARRAY<T>& energy_gradient,const ARRAY<T>& invariants,const int id) const=0;
    virtual void Energy_Hessian(ARRAY<T>& energy_hessian,const ARRAY<T>& invariants,const int id) const=0;
//#####################################################################
};
}
#endif
