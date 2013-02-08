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

namespace PhysBAM{

using ::std::exp

template<class TV>
class MPM_CONSTITUTIVE_MODEL
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef MATRIX<T,TV::m> MTX;
    typedef DIAGONAL_MATRIX<T,TV::m> MTX_DIAG;

public:

    T mu0,lambda0,xi;

    MTX Fe,Fp;

    T Je,Jp;
    MTX Re,Se,Ue,Ve;
    MTX_DIAG SIGMAe;
    T mu,lambda;

    MPM_CONSTITUTIVE_MODEL();
    ~MPM_CONSTITUTIVE_MODEL();
    void Initialize(const T youngs_modulus,const T poisson_ratio,const T hardening_coefficient);

    void Update_Quantities_Using_Current_Deformation_Gradient();

    T Energy_Density_Psi();
    
    MTX dPsi_dFe();

    MTX d2Psi_dFe_dFe_Action_dF(const MTX& dF);

//#####################################################################
};
}
#endif
