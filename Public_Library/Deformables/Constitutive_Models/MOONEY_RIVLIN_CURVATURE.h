//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
// Class MOONEY_RIVLIN_CURVATURE
//#####################################################################
#ifndef __MOONEY_RIVLIN_CURVATURE__
#define __MOONEY_RIVLIN_CURVATURE__
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>

namespace PhysBAM{

template<class T>
class MOONEY_RIVLIN_CURVATURE
{
    typedef VECTOR<T,3> TV;
    typedef MATRIX<T,3> TM;
    typedef SYMMETRIC_MATRIX<T,3> SM;
    typedef MATRIX<T,3,5> TM2;
    typedef MATRIX<MATRIX<T,3,3>,5,5> TT;

public:
    const T c1;
    const T c2;
    const T thickness;

    MOONEY_RIVLIN_CURVATURE(T c1,T c2,T thickness);

    //#####################################################################
    T Potential_Energy(TV A1,TV A2,TV dA11,TV dA12,TV dA22,const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge, TT& he) const;
    T Potential_Energy(TV A1,TV A2,TV dA11,TV dA12,TV dA22,const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge) const;
};
}
#endif
