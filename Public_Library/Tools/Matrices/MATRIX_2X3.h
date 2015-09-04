//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_2X3
//#####################################################################
#ifndef __MATRIX_2X3__
#define __MATRIX_2X3__

#include <Tools/Matrices/MATRIX_BASE.h>
namespace PhysBAM{

template<class T_input>
class MATRIX<T_input,2,3>:public MATRIX_BASE<T_input,MATRIX<T_input,2,3> >
{
public:
    typedef T_input T;typedef T SCALAR;
    typedef MATRIX_BASE<T,MATRIX<T,2,3> > BASE;
    enum WORKAROUND1 {m=2,n=3};
    using BASE::operator*;

    MATRIX<T,3,2> transpose;

    MATRIX()
    {}

    MATRIX(INITIAL_SIZE m,INITIAL_SIZE n)
    {
        assert(m==INITIAL_SIZE(2) && n==INITIAL_SIZE(3));
    }

    MATRIX(const MATRIX& matrix)
        :transpose(matrix.transpose)
    {}

    template<class T2>
    MATRIX(const MATRIX<T2,2,3>& matrix)
        :transpose(MATRIX<T,3,2>(matrix.transpose))
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

    MATRIX<T,2,3>& operator=(const MATRIX<T,2,3>& M)
    {
        transpose=M.transpose;
        return *this;
    }

    int Rows() const
    {return transpose.Columns();}

    int Columns() const
    {return transpose.Rows();}

    T& operator()(const int i,const int j)
    {return transpose(j,i);}

    const T& operator()(const int i,const int j) const
    {return transpose(j,i);}

    bool Valid_Index(const int i,const int j) const
    {return transpose.Valid_Index(j,i);}

    bool operator==(const MATRIX<T,2,3>& A) const
    {return transpose==A.transpose;}

    bool operator!=(const MATRIX<T,2,3>& A) const
    {return transpose!=A.transpose;}

    MATRIX<T,2,3> operator-() const
    {MATRIX<T,2,3> r;r.transpose=-transpose;return r;}

    MATRIX<T,2,3> operator*(const T a) const
    {MATRIX<T,2,3> r;r.transpose=transpose*a;return r;}

    MATRIX<T,2,3>& operator*=(const T a)
    {transpose*=a;return *this;}

    MATRIX<T,2,3> operator/(const T a) const
    {MATRIX<T,2,3> r;r.transpose=transpose/a;return r;}

    MATRIX<T,2,3>& operator/=(const T a)
    {transpose/=a;return *this;}

    MATRIX<T,2,3> operator+(const MATRIX<T,2,3>& A) const
    {MATRIX<T,2,3> r;r.transpose=transpose+A.transpose;return r;}

    MATRIX<T,2,3> operator+=(const MATRIX<T,2,3>& A)
    {transpose+=A.transpose;return *this;}

    MATRIX<T,2,3> operator-(const MATRIX<T,2,3>& A) const
    {MATRIX<T,2,3> r;r.transpose=transpose-A.transpose;return r;}

    MATRIX<T,2,3> operator-=(const MATRIX<T,2,3>& A)
    {transpose-=A.transpose;return *this;}

    // template<class TM>
    // auto operator*(const TM& A) const
    // {return transpose.Transpose_Times(A);}

    static MATRIX<T,2,3> Transposed(const MATRIX<T,3,2>& A)
    {MATRIX<T,2,3> transpose;transpose.transpose=A;return transpose;} 

    const MATRIX<T,3,2>& Transposed() const
    {return transpose;}

    MATRIX<T,2,3> Cofactor_Matrix() const
    {return Transposed(transpose.Cofactor_Matrix());}

    static MATRIX<T,2,3> Outer_Product(const VECTOR<T,m>& u,const VECTOR<T,n>& v)
    {MATRIX<T,2,3> result;result.transpose=MATRIX<T,3,2>::Outer_Product(v,u);return result;}

    T Max_Abs() const
    {return transpose.Max_Abs();}

    T Frobenius_Norm_Squared() const
    {return transpose.Frobenius_Norm_Squared();}

    template<class TM> auto Transpose_Times(const TM& A) const
    {return transpose*A;}

    template<class TM> auto Times_Transpose(const TM& A) const
    {return (A*transpose).Transposed();}

    VECTOR<T,2> Column(const int j) const
    {return transpose.Row(j);}

    void Set_Column(const int j,const VECTOR<T,2>& v)
    {transpose.Set_Row(j,v);}

    void Add_Column(const int j,const VECTOR<T,2>& v)
    {transpose.Add_Row(j,v);}

    VECTOR<T,3> Row(const int j) const
    {return transpose.Column(j);}

    void Set_Row(const int j,const VECTOR<T,3>& v)
    {transpose.Set_Column(j,v);}

    void Add_Row(const int j,const VECTOR<T,3>& v)
    {transpose.Add_Column(j,v);}

    void Fast_Singular_Value_Decomposition(MATRIX<T,2>& U,DIAGONAL_MATRIX<T,2>& singular_values,MATRIX<T,3,2>& V) const
    {transpose.Fast_Singular_Value_Decomposition(V,singular_values,U);}

//#####################################################################
};

template<class T> MATRIX<T,2,3>
operator*(const T& a,const MATRIX<T,2,3>& A)
{
    MATRIX<T,2,3> result;
    result.transpose=A.transpose*a;
    return result;
}

// template<class TM0,class T>
// typename enable_if<IS_MATRIX<TM0>::value,decltype(TM0().Times_Transpose(TM1()))>::type
// operator*(const TM0& A,const MATRIX<T,2,3>& B)
// {
//     return A.Times_Transpose(B.transpose);
// }
}
#endif
