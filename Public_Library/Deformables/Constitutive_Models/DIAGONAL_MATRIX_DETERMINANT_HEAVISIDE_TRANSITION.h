//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONAL_MATRIX_DETERMINANT_HEAVISIDE_TRANSITION
//#####################################################################
#ifndef __DIAGONAL_MATRIX_DETERMINANT_HEAVISIDE_TRANSITION__
#define __DIAGONAL_MATRIX_DETERMINANT_HEAVISIDE_TRANSITION__

#include <Tools/Log/DEBUG_UTILITIES.h>

namespace PhysBAM{

template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d> class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE;
template<class T,int d, class TRANSITION_T>
class DIAGONAL_MATRIX_DETERMINANT_HEAVISIDE_TRANSITION
{

private:

    TRANSITION_T base;

    void DDH_Helper(const DIAGONAL_MATRIX<T,2>& S,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& d2HdS2) const
    {
        T J   = S.Determinant();
        T HJ  = base.Hx(J);
        T HJJ = base.Hxx(J);
        
        d2HdS2.x1111 = HJJ*sqr(S.x.y);
        d2HdS2.x2222 = HJJ*sqr(S.x.x);
        d2HdS2.x2211 = HJJ*S.x.x*S.x.y + HJ;
        
        d2HdS2.x2112 = -HJ;
        d2HdS2.x2121 = 0;
    }

    void DDH_Helper(const DIAGONAL_MATRIX<T,3>& S,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& d2HdS2) const
    {
        T J   = S.Determinant();
        T HJ  = base.Hx(J);
        T HJJ = base.Hxx(J);
        
        d2HdS2.x1111 = HJJ*sqr(S.x.y*S.x.z);
        d2HdS2.x2222 = HJJ*sqr(S.x.x*S.x.z);
        d2HdS2.x3333 = HJJ*sqr(S.x.x*S.x.y);
 
        d2HdS2.x2211 = HJJ*S.x.x*sqr(S.x.z)*S.x.y + HJ*S.x.z;
        d2HdS2.x3311 = HJJ*S.x.x*sqr(S.x.y)*S.x.z + HJ*S.x.y;
        d2HdS2.x3322 = HJJ*S.x.z*sqr(S.x.x)*S.x.y + HJ*S.x.x;
        
        d2HdS2.x2112 = -HJ*S.x.z;
        d2HdS2.x3113 = -HJ*S.x.y;
        d2HdS2.x3223 = -HJ*S.x.x;

        d2HdS2.x2121 = 0; 
        d2HdS2.x3131 = 0; 
        d2HdS2.x3232 = 0; 
    }

public:

    DIAGONAL_MATRIX_DETERMINANT_HEAVISIDE_TRANSITION() {}

    DIAGONAL_MATRIX_DETERMINANT_HEAVISIDE_TRANSITION(const T input_J_min, const T input_J_max)
        :base(input_J_min,input_J_max)
    {}

    ~DIAGONAL_MATRIX_DETERMINANT_HEAVISIDE_TRANSITION() {}
    
    void Initialize(const T input_J_min, const T input_J_max)
    {return base.Initialize(input_J_min,input_J_max);}

    T H(const DIAGONAL_MATRIX<T,d>& S) const
    {return base.H(S.Determinant());}

    DIAGONAL_MATRIX<T,d> DH(const DIAGONAL_MATRIX<T,d>& S) const
    {return base.Hx(S.Determinant())*S.Cofactor_Matrix();}

    void DDH(const DIAGONAL_MATRIX<T,d>& S,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& d2HdS2) const
    {return DDH_Helper(S,d2HdS2);}
};
}

#endif
