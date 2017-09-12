//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE
//#####################################################################
#ifndef __DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE__
#define __DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE__

#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

namespace{
template<class T> MATRIX<T,3>
Differential_Helper(const SYMMETRIC_MATRIX<T,3>& H,const VECTOR<T,3>& B,const VECTOR<T,3>& C,const MATRIX<T,3>& dF)
{
    return MATRIX<T,3>(H(0,0)*dF.x[0]+H(1,0)*dF.x[4]+H(2,0)*dF.x[8],B(2)*dF.x[1]+C(2)*dF.x[3],B(1)*dF.x[2]+C(1)*dF.x[6],
        C(2)*dF.x[1]+B(2)*dF.x[3],H(1,0)*dF.x[0]+H(1,1)*dF.x[4]+H(2,1)*dF.x[8],B(0)*dF.x[5]+C(0)*dF.x[7],
        C(1)*dF.x[2]+B(1)*dF.x[6],C(0)*dF.x[5]+B(0)*dF.x[7],H(2,0)*dF.x[0]+H(2,1)*dF.x[4]+H(2,2)*dF.x[8]);
}
template<class T> MATRIX<T,2>
Differential_Helper(const SYMMETRIC_MATRIX<T,2>& H,const VECTOR<T,1>& B,const VECTOR<T,1>& C,const MATRIX<T,2>& dF)
{
    return MATRIX<T,2>(H(0,0)*dF.x[0]+H(1,0)*dF.x[3],B(10)*dF.x[1]+C(10)*dF.x[2],C(10)*dF.x[1]+B(10)*dF.x[2],H(1,0)*dF.x[0]+H(1,1)*dF.x[3]);
}
template<class T> MATRIX<T,1>
Differential_Helper(const SYMMETRIC_MATRIX<T,1>& H,const VECTOR<T,0>& B,const VECTOR<T,0>& C,const MATRIX<T,1>& dF)
{
    return MATRIX<T,1>(H(0,0)*dF.x[0]);
}
}

template<class TV>
class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE
{
    typedef typename TV::SCALAR T;
public:
    SYMMETRIC_MATRIX<T,TV::m> H; // H(i,k) = x_iikk
    typename TV::SPIN B,C; // B = x_ikik; C = x_ikki; order: (1D: none; 2D: 01; 3D: 12 20 01)

    MATRIX<T,TV::m> Differential(const MATRIX<T,TV::m>& dF) const
    {return Differential_Helper(H,B,C,dF);}

//#####################################################################
    void Enforce_Definiteness();
    void Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,TV::m>& F,const TV& dE_ds,const SYMMETRIC_MATRIX<T,TV::m>& dE_dsds); // Not robust
//#####################################################################
};
}
#endif
