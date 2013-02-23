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

template<class TV>
class MPM_CONSTITUTIVE_MODEL
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_CONSTITUTIVE_MODEL();
    ~MPM_CONSTITUTIVE_MODEL();

    void Compute_Helper_Quantities_Using_F(const MATRIX<T,TV::m>& Fe,const MATRIX<T,TV::m>& Fp,T& Je,MATRIX<T,TV::m>& Re,MATRIX<T,TV::m>& Se) const;
    T Compute_Elastic_Energy_Density_Psi(const T& mu,const T& lambda,const MATRIX<T,TV::m>& Fe,const MATRIX<T,TV::m>& Re,const T& Je) const;
    MATRIX<T,TV::m> Compute_dPsi_dFe(const T& mu,const T& lambda,const MATRIX<T,TV::m>& Fe,const MATRIX<T,TV::m>& Re,const T& Je) const;
    MATRIX<T,TV::m> Compute_d2Psi_dFe_dFe_Action_dF(const T& mu,const T& lambda,const MATRIX<T,TV::m>& Fe,const T& Je,const MATRIX<T,TV::m>& Re,const MATRIX<T,TV::m>& Se,const MATRIX<T,TV::m>& dF) const;
    void Derivative_Test() const;
//#####################################################################
};
}
#endif
