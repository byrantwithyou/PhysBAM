//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_2X3
//#####################################################################
#ifndef __MATRIX_2X3__
#define __MATRIX_2X3__

#include <Tools/Matrices/TRANSPOSE_MATRIX.h>
namespace PhysBAM{

template<class T_input>
class MATRIX<T_input,2,3>:public TRANSPOSE_MATRIX<MATRIX<T_input,3,2> >
{
public:
    typedef T_input T;typedef T SCALAR;
    typedef TRANSPOSE_MATRIX<MATRIX<T_input,3,2> > BASE;
    using BASE::transpose;

    MATRIX()
    {}

    MATRIX(INITIAL_SIZE m,INITIAL_SIZE n)
    {
        assert(m==INITIAL_SIZE(2) && n==INITIAL_SIZE(3));
    }

    MATRIX(const MATRIX& matrix)
        :TRANSPOSE_MATRIX<MATRIX<T,3,2> >(static_cast<const TRANSPOSE_MATRIX<MATRIX<T,3,2> >&>(matrix))
    {}

    template<class T2>
    MATRIX(const MATRIX<T2,2,3>& matrix)
        :TRANSPOSE_MATRIX<MATRIX<T,3,2> >(static_cast<const TRANSPOSE_MATRIX<MATRIX<T2,3,2> >&>(matrix))
    {}

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        assert(A.Rows()==2 && A.Columns()==3);for(int j=0;j<3;j++) for(int i=0;i<2;i++) (*this)(i,j)=A(i,j);
    }

    MATRIX(const T x00,const T x10,const T x01,const T x11,const T x02,const T x12)
    {
        transpose=MATRIX<T,3,2>(x00,x01,x02,x10,x11,x12);
    }

    MATRIX(const VECTOR<T,2>& v1,const VECTOR<T,2>& v2,const VECTOR<T,2>& v3)
    {
        transpose=MATRIX<T,3,2>(v1.x,v2.x,v3.x,v1.y,v2.y,v3.y);
    }

    void Fast_Singular_Value_Decomposition(MATRIX<T,2>& U,DIAGONAL_MATRIX<T,2>& singular_values,MATRIX<T,3,2>& V) const
    {transpose.Fast_Singular_Value_Decomposition(V,singular_values,U);}

//#####################################################################
};
}
#endif
