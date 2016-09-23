//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
// Class MOONEY_RIVLIN_CURVATURE
//#####################################################################
#ifndef __MOONEY_RIVLIN_CURVATURE__
#define __MOONEY_RIVLIN_CURVATURE__
#include <Core/Matrices/SYMMETRIC_MATRIX.h>

namespace PhysBAM{

template<class T>
class MOONEY_RIVLIN_CURVATURE
{
    typedef VECTOR<T,3> TV;
    typedef MATRIX<T,3> TM;
    typedef SYMMETRIC_MATRIX<T,3> SM;
    typedef VECTOR<VECTOR<T,3>,5> TM2;
    typedef MATRIX<MATRIX<T,3,3>,5,5> TT;

public:
    const T c1;
    const T c2;
    const T thickness;

    MOONEY_RIVLIN_CURVATURE(T c1,T c2,T thickness);

    //#####################################################################
    T Potential_Energy(TV A1,TV A2,TV dA11,TV dA12,TV dA22,const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge, TT& he,T weight=1) const;
    T Potential_Energy(TV A1,TV A2,TV dA11,TV dA12,TV dA22,const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge,T weight=1) const;
    template<class T1,class T2,class T3,class T4,class T5,class T6> T
    Helper(const T1& a1,const T2& a2,const T3& lambda_n,const T4& dlambda_n1,const T5& dlambda_n2,
        const TM& S,TM2& ge,T6 he,T w,T scale) const;
    template<class T1,class T2,class T3,class T4,class T5,class T6> T
    Potential_Energy_Helper(const T1& a1,const T2& a2,const T3& da11,const T4& da12,const T5& da22,
        const VECTOR<TM,3>& G0_inv,const VECTOR<T,3>& G0_det,TM2& ge,T6 he,T weight) const;
    template<class T1> void Extract_Helper(const T1& in,TT* he) const;
    template<class T1> void Extract_Helper(const T1& in,int he) const;
};
}
#endif
