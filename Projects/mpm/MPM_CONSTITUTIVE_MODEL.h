//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_CONSTITUTIVE_MODEL
//#####################################################################
#ifndef __MPM_CONSTITUTIVE_MODEL__
#define __MPM_CONSTITUTIVE_MODEL__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>

namespace PhysBAM{

using ::std::exp;

template<class TV>
class MPM_CONSTITUTIVE_MODEL
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:

    T mu0,lambda0,xi;

    MATRIX<T,TV::m> Fe,Fp;

    T Je,Jp;
    MATRIX<T,TV::m> Re,Se,Ue,Ve;
    DIAGONAL_MATRIX<T,TV::m> SIGMAe;
    T mu,lambda;

    MPM_CONSTITUTIVE_MODEL();
    ~MPM_CONSTITUTIVE_MODEL();
    void Initialize(const T youngs_modulus,const T poisson_ratio,const T hardening_coefficient);

    void Update_Quantities_Using_Current_Deformation_Gradient();

    T Energy_Density_Psi();
    
    MATRIX<T,TV::m> dPsi_dFe();

    MATRIX<T,TV::m> d2Psi_dFe_dFe_Action_dF(const MATRIX<T,TV::m>& dF);

protected:
//#####################################################################
};
}
#endif
