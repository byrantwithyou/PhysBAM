//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_VIEW
//#####################################################################
#ifndef __MATRIX_VIEW__
#define __MATRIX_VIEW__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/MATRIX_BASE.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <type_traits>
namespace PhysBAM{

template<class T>
class MATRIX_VIEW:public MATRIX_BASE<typename remove_const<T>::type,MATRIX_VIEW<T> >
{
    template<class S> struct COPY_CONST{typedef typename conditional<is_const<T>::value,typename add_const<S>::type,S>::type TYPE;};
    struct UNUSABLE{};
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef MATRIX_BASE<T,MATRIX_VIEW<T> > BASE;
    typedef typename remove_const<T>::type SCALAR;
    typedef typename conditional<is_const<T>::value,SCALAR,UNUSABLE>::type CT;
    typedef typename conditional<is_const<T>::value,UNUSABLE,SCALAR>::type NT;
    int m,n,s; // size of the m by n matrix, s is stride
    T* x; // one dimensional data

    MATRIX_VIEW()
        :m(0),n(0),s(0),x(0)
    {}

    MATRIX_VIEW(T* x,int m,int n)
        :m(m),n(n),s(m),x(m*n)
    {
    }

    MATRIX_VIEW(T* x,int m,int n,int s)
        :m(m),n(n),s(s),x(m*n)
    {
    }

    template<int a,int b>
    MATRIX_VIEW(const MATRIX<CT,a,b>& M)
        :m(a),n(b),s(a),x(M.x)
    {
    }

    template<int a,int b>
    MATRIX_VIEW(MATRIX<NT,a,b>& M)
        :m(a),n(b),s(a),x(M.x)
    {
    }

    MATRIX_VIEW(const MATRIX_MXN<CT>& M)
        :m(M.m),n(M.n),s(M.m),x(M.x.base_pointer)
    {
    }

    MATRIX_VIEW(MATRIX_MXN<NT>& M)
        :m(M.m),n(M.n),s(M.m),x(M.x.base_pointer)
    {
    }

    MATRIX_VIEW(const MATRIX_VIEW<T>& A) = default;
    MATRIX_VIEW(MATRIX_VIEW&& array) = default;
    ~MATRIX_VIEW() = default;

    template<class T_MATRIX>
    MATRIX_VIEW& operator=(const MATRIX_BASE<T,T_MATRIX>& A)
    {
        if((void*)&A==(void*)this) return *this;
        assert(m==A.Rows());
        assert(n==A.Columns());
        for(int i=0;i<m;i++)
            for(int j=0;j<n;j++)
                (*this)(i,j)=A(i,j);
        return *this;
    }

    MATRIX_VIEW& operator=(const MATRIX_VIEW& source)
    {return *this=static_cast<const BASE&>(source);}

    MATRIX_VIEW& operator=(MATRIX_VIEW&& source) = delete;

    void Set(T* raw_data,int m_input,int n_input)
    {m=m_input;n=n_input;s=m_input;x=raw_data;}

    void Set(T* raw_data,int m_input,int n_input,int s_input)
    {m=m_input;n=n_input;s=s_input;x=raw_data;}

    void Set(MATRIX_VIEW<T> array)
    {Set(array.x,array.m,array.n,array.s);}

    template<int a,int b>
    void Set(MATRIX<NT,a,b>& array)
    {Set(array.x,a,b);}

    template<int a,int b>
    void Set(const MATRIX<CT,a,b>& array)
    {Set(array.x,a,b);}

    void Set(MATRIX_MXN<NT>& array)
    {Set(array.x.Get_Array_Pointer(),array.m,array.n);}

    void Set(const MATRIX_MXN<CT>& array)
    {Set(array.x.Get_Array_Pointer(),array.m,array.n);}

    int Rows() const
    {return m;}

    int Columns() const
    {return n;}

    T& operator()(const int i,const int j)
    {assert((unsigned)i<(unsigned)m);assert((unsigned)j<(unsigned)n);return x[j*s+i];}

    const T& operator()(const int i,const int j) const
    {assert((unsigned)i<(unsigned)m);assert((unsigned)j<(unsigned)n);return x[j*s+i];}

    bool Valid_Index(const int i,const int j) const
    {return (unsigned)i<(unsigned)m && (unsigned)j<(unsigned)n;}

//#####################################################################
};
}
#endif
