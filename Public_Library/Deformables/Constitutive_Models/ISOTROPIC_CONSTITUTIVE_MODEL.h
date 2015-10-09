//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ISOTROPIC_CONSTITUTIVE_MODEL
//##################################################################### 
#ifndef __ISOTROPIC_CONSTITUTIVE_MODEL__
#define __ISOTROPIC_CONSTITUTIVE_MODEL__

#include <Deformables/Constitutive_Models/CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T,int d> class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE;

template<class T,int d>
class ISOTROPIC_CONSTITUTIVE_MODEL:public CONSTITUTIVE_MODEL<T,d>
{
public:
    using CONSTITUTIVE_MODEL<T,d>::Isotropic_Stress_Derivative;
    ISOTROPIC_CONSTITUTIVE_MODEL();
    virtual ~ISOTROPIC_CONSTITUTIVE_MODEL();

//#####################################################################
    virtual MATRIX<T,d> dP_From_dF(const MATRIX<T,d>& dF,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dPi_dF,const int id) const;
    virtual DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const=0;
    virtual void Update_State_Dependent_Auxiliary_Variables(const DIAGONAL_MATRIX<T,d>& F,const int id);
    virtual T Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const;
    void Test(const DIAGONAL_MATRIX<T,d>& F,const int id) const;
//#####################################################################
};
}
#endif
