//#####################################################################
// Copyright 2007, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_BASE
//#####################################################################
#ifndef __MATRIX_BASE__
#define __MATRIX_BASE__

#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Data_Structures/ELEMENT_ID.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/max.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/MATRIX_FORWARD.h>
#include <cmath>
#include <iomanip>
#include <iostream>
using std::abs;
namespace PhysBAM{

template<class T_MATRIX> struct MATRIX_INFO;
template<class T,int a,int b> struct MATRIX_INFO<MATRIX<T,a,b> >{static const int m=a,n=b;};
template<class T,int d> struct MATRIX_INFO<DIAGONAL_MATRIX<T,d> >{static const int m=d,n=d;};
template<class T,int d> struct MATRIX_INFO<SYMMETRIC_MATRIX<T,d> >{static const int m=d,n=d;};
template<class T,int d> struct MATRIX_INFO<UPPER_TRIANGULAR_MATRIX<T,d> >{static const int m=d,n=d;};
template<class T> struct MATRIX_INFO<MATRIX_MXN<T> >{static const int m=-1,n=-1;};
template<class T> struct MATRIX_INFO<MATRIX_VIEW<T> >{static const int m=-1,n=-1;};
template<class T,class T_MATRIX> struct MATRIX_INFO<MATRIX_BASE<T,T_MATRIX> >:public MATRIX_INFO<T_MATRIX> {};

template<class T,int a,int b> struct SIZED_MATRIX_HELPER {typedef MATRIX<T,a,b> TYPE;};
template<class T,int a> struct SIZED_MATRIX_HELPER<T,a,-1> {typedef MATRIX_MXN<T> TYPE;};
template<class T,int b> struct SIZED_MATRIX_HELPER<T,-1,b> {typedef MATRIX_MXN<T> TYPE;};
template<class T> struct SIZED_MATRIX_HELPER<T,-1,-1> {typedef MATRIX_MXN<T> TYPE;};
template<class T,int a,int b> using SIZED_MATRIX = typename SIZED_MATRIX_HELPER<T,a,b>::TYPE;
template<class T_MAT> using MATRIX_COPY = SIZED_MATRIX<typename T_MAT::ELEMENT,MATRIX_INFO<T_MAT>::m,MATRIX_INFO<T_MAT>::n>;

template<class T,int a> struct SIZED_VECTOR_HELPER {typedef VECTOR<T,a> TYPE;};
template<class T> struct SIZED_VECTOR_HELPER<T,-1> {typedef ARRAY<T> TYPE;};
template<class T,int a> using SIZED_VECTOR = typename SIZED_VECTOR_HELPER<T,a>::TYPE;
template<class T_VEC> using VECTOR_COPY = SIZED_VECTOR<typename T_VEC::ELEMENT,FIXED_SIZE_VECTOR<T_VEC>::size>;

template<class T_MATRIX> struct EFFICIENT_MATRIX {static const bool value=false;};
template<class T,int m,int n> struct EFFICIENT_MATRIX<MATRIX<T,m,n> > {static const bool value=((m>=2 && m<=3 && n>=2 && n<=3) || (m==4 && n==4) || (m==0 && n==0));};
template<class T,int d> struct EFFICIENT_MATRIX<DIAGONAL_MATRIX<T,d> > {static const bool value=(d==2 || d==3);};
template<class T,int d> struct EFFICIENT_MATRIX<SYMMETRIC_MATRIX<T,d> > {static const bool value=(d==2 || d==3);};
template<class T,class T_MATRIX> struct EFFICIENT_MATRIX<MATRIX_BASE<T,T_MATRIX> >:public EFFICIENT_MATRIX<T_MATRIX>{};
template<class T,int d> struct EFFICIENT_MATRIX<VECTOR<T,d> > {static const bool value=(d<=3);};
template<class T,class T_VECTOR> struct EFFICIENT_MATRIX<ARRAY_BASE<T,T_VECTOR> >:public EFFICIENT_MATRIX<T_VECTOR>{};


namespace{
template<int line,class A,class B=void> struct ASSERT_EFFICIENT
{
    template<class T> struct EFFICIENT_OR_VOID {static const bool value=EFFICIENT_MATRIX<T>::value || is_same<T,void>::value;};
    static const bool efficient=EFFICIENT_OR_VOID<A>::value && EFFICIENT_OR_VOID<B>::value;
    struct UNUSABLE{};

    ASSERT_EFFICIENT(typename conditional<efficient,UNUSABLE,const char*>::type str)
    {}

    ASSERT_EFFICIENT(typename conditional<efficient,const char*,UNUSABLE>::type str)
    {PHYSBAM_WARNING(std::string("Base implementation used for: ")+str+".");}
};

#ifdef NDEBUG
#define WARN_IF_NOT_EFFICIENT(...)
#else
#ifdef _WIN32
#define WARN_IF_NOT_EFFICIENT(...) ASSERT_EFFICIENT<__LINE__,__VA_ARGS__>(__FUNCTION__)
#else
#define WARN_IF_NOT_EFFICIENT(...) ASSERT_EFFICIENT<__LINE__,__VA_ARGS__>(__PRETTY_FUNCTION__)
#endif
#endif
}

template<class T,class T_MATRIX>
class MATRIX_BASE
{
protected:
    MATRIX_BASE& operator=(const MATRIX_BASE&) = default;
    MATRIX_BASE& operator=(MATRIX_BASE&&) = default;
    MATRIX_BASE(const MATRIX_BASE&) = default;

public:
    typedef T SCALAR;
    typedef SIZED_VECTOR<T,MATRIX_INFO<T_MATRIX>::m> LEFT_VECTOR;
    typedef SIZED_VECTOR<T,MATRIX_INFO<T_MATRIX>::n> RIGHT_VECTOR;
    typedef SIZED_VECTOR<int,MATRIX_INFO<T_MATRIX>::n> COLUMN_PERMUTATION;

    MATRIX_BASE() = default;

    ~MATRIX_BASE() = default;

    T_MATRIX& Derived()
    {return static_cast<T_MATRIX&>(*this);}

    const T_MATRIX& Derived() const
    {return static_cast<const T_MATRIX&>(*this);}

    T operator()(const int i,const int j) const
    {return Derived()(i,j);}

    T& operator()(const int i,const int j)
    {return Derived()(i,j);}

    int Rows() const
    {return Derived().Rows();}

    int Columns() const
    {return Derived().Columns();}

protected:
    // Logic only; no type or bounds checking; requires all arguments to be distinct
    template<class T_MATRIX0,class T_MATRIX1,class T_MATRIX2>
    static void Add_Times_Matrix_Helper(const T_MATRIX0& A,const T_MATRIX1& B,T_MATRIX2& C)
    {for(int j=0;j<B.Columns();j++) for(int k=0;k<A.Columns();k++) for(int i=0;i<A.Rows();i++) C(i,j)+=A(i,k)*B(k,j);}

    template<class T_MATRIX0,class T_MATRIX1,class T_MATRIX2>
    static void Add_Transpose_Times_Matrix_Helper(const T_MATRIX0& A,const T_MATRIX1& B,T_MATRIX2& C)
    {for(int j=0;j<B.Columns();j++) for(int i=0;i<A.Columns();i++) for(int k=0;k<A.Rows();k++) C(i,j)+=A(k,i)*B(k,j);}

    template<class T_MATRIX0,class T_MATRIX1,class T_MATRIX2>
    static void Add_Times_Transpose_Matrix_Helper(const T_MATRIX0& A,const T_MATRIX1& B,T_MATRIX2& C)
    {for(int k=0;k<A.Columns();k++) for(int j=0;j<B.Rows();j++) for(int i=0;i<A.Rows();i++) C(i,j)+=A(i,k)*B(j,k);}

    template<class T_MATRIX0,class T_VECTOR,class T_VECTOR1>
    static void Add_Transpose_Times_Vector_Helper(const T_MATRIX0& A,const ARRAY_BASE<T,T_VECTOR>& v,ARRAY_BASE<T,T_VECTOR1>& y)
    {for(int i=0;i<y.Size();i++) for(int j=0;j<v.Size();j++) y(i)+=A(j,i)*v(j);}

    template<class T_MATRIX0,class T_VECTOR,class T_VECTOR1>
    static void Add_Times_Vector_Helper(const T_MATRIX0& A,const ARRAY_BASE<T,T_VECTOR>& v,ARRAY_BASE<T,T_VECTOR1>& y)
    {for(int j=0;j<v.Size();j++) for(int i=0;i<y.Size();i++) y(i)+=A(i,j)*v(j);}

    template<class T_MATRIX0,class T_VECTOR,class T_VECTOR1>
    static void Subtract_Times_Vector_Helper(const T_MATRIX0& A,const ARRAY_BASE<T,T_VECTOR>& v,ARRAY_BASE<T,T_VECTOR1>& y)
    {for(int j=0;j<v.Size();j++) for(int i=0;i<y.Size();i++) y(i)-=A(i,j)*v(j);}

    template<class T_MATRIX0,class T_MATRIX1>
    static bool Test_Aliased(const T_MATRIX0& A,const T_MATRIX1& B)
    {assert((void*)&A!=(void*)&B);return false;}

    template<class T_MATRIX0>
    static bool Test_Aliased(const T_MATRIX0& A,const T_MATRIX0& B)
    {return &A==&B;}
public:

    // With bounds checking; arguments may be the same.
    template<class T_MATRIX0,class T_MATRIX1,class T_MATRIX2>
    static void Add_Times(const T_MATRIX0& A,const T_MATRIX1& B,MATRIX_BASE<T,T_MATRIX2>& C)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX0,T_MATRIX1);assert(A.Columns()==B.Rows() && C.Rows()==A.Rows() && C.Columns()==B.Columns());
    if(Test_Aliased(A,C) || Test_Aliased(B,C)){T_MATRIX2 M(INITIAL_SIZE(C.Rows()),INITIAL_SIZE(C.Columns()));
    Add_Times_Matrix_Helper(A,B,M);C+=M;}else Add_Times_Matrix_Helper(A,B,C);}

    template<class T_MATRIX0,class T_MATRIX1,class T_MATRIX2>
    static void Add_Transpose_Times(const T_MATRIX0& A,const T_MATRIX1& B,MATRIX_BASE<T,T_MATRIX2>& C)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX0,T_MATRIX1);assert(A.Rows()==B.Rows() && C.Rows()==A.Columns() && C.Columns()==B.Columns());
    if(Test_Aliased(A,C) || Test_Aliased(B,C)){T_MATRIX2 M(INITIAL_SIZE(C.Rows()),INITIAL_SIZE(C.Columns()));
    Add_Transpose_Times_Matrix_Helper(A,B,M);C+=M;}else Add_Transpose_Times_Matrix_Helper(A,B,C);}

    template<class T_MATRIX0,class T_MATRIX1,class T_MATRIX2>
    static void Add_Times_Transpose(const T_MATRIX0& A,const T_MATRIX1& B,MATRIX_BASE<T,T_MATRIX2>& C)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX0,T_MATRIX1);assert(A.Columns()==B.Columns() && C.Rows()==A.Rows() && C.Columns()==B.Rows());
    if(Test_Aliased(A,C) || Test_Aliased(B,C)){T_MATRIX2 M(INITIAL_SIZE(C.Rows()),INITIAL_SIZE(C.Columns()));
    Add_Times_Transpose_Matrix_Helper(A,B,M);C+=M;}else Add_Times_Transpose_Matrix_Helper(A,B,C);}

    template<class T_MATRIX0,class T_VECTOR,class T_VECTOR1>
    static void Add_Transpose_Times(const MATRIX_BASE<T,T_MATRIX0>& A,const ARRAY_BASE<T,T_VECTOR>& v,ARRAY_BASE<T,T_VECTOR1>& y)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX0,T_VECTOR1);assert(A.Rows()==v.Size() && A.Columns()==y.Size());
    if(Test_Aliased(v,y)){VECTOR_COPY<T_VECTOR> u(v);Add_Transpose_Times_Vector_Helper(A,u,y);}else Add_Transpose_Times_Vector_Helper(A,v,y);}

    template<class T_MATRIX0,class T_VECTOR,class T_VECTOR1>
    static void Add_Times(const MATRIX_BASE<T,T_MATRIX0>& A,const ARRAY_BASE<T,T_VECTOR>& v,ARRAY_BASE<T,T_VECTOR1>& y)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX0,T_VECTOR);assert(A.Columns()==v.Size() && A.Rows()==y.Size());
    if(Test_Aliased(v,y)){VECTOR_COPY<T_VECTOR> u(v);Add_Times_Vector_Helper(A,u,y);}else Add_Times_Vector_Helper(A,v,y);}

    template<class T_MATRIX0,class T_VECTOR,class T_VECTOR1>
    static void Subtract_Times(const MATRIX_BASE<T,T_MATRIX0>& A,const ARRAY_BASE<T,T_VECTOR>& v,ARRAY_BASE<T,T_VECTOR1>& y)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX0,T_VECTOR);assert(A.Columns()==v.Size() && A.Rows()==y.Size());
    if(Test_Aliased(v,y)){VECTOR_COPY<T_VECTOR> u(v);Subtract_Times_Vector_Helper(A,u,y);}else Subtract_Times_Vector_Helper(A,v,y);}

    template<class T_MATRIX0>
    auto Transpose_Times(const MATRIX_BASE<T,T_MATRIX0>& A) const
    {
        static const int a=MATRIX_INFO<T_MATRIX>::n;
        static const int b=MATRIX_INFO<T_MATRIX0>::n;
        SIZED_MATRIX<T,a,b> M((INITIAL_SIZE)Columns(),(INITIAL_SIZE)A.Columns());
        Add_Transpose_Times(Derived(),A,M);
        return M;
    }

    template<int d>
    auto Transpose_Times(const SYMMETRIC_MATRIX<T,d>& A) const
    {
        static const int a=MATRIX_INFO<T_MATRIX>::n;
        SIZED_MATRIX<T,a,d> M((INITIAL_SIZE)Columns(),(INITIAL_SIZE)d);
        Add_Transpose_Times(Derived(),A,M);
        return M;
    }

    template<int d>
    auto Transpose_Times(const DIAGONAL_MATRIX<T,d>& A) const
    {
        WARN_IF_NOT_EFFICIENT(T_MATRIX);
        static const int a=MATRIX_INFO<T_MATRIX>::n;
        SIZED_MATRIX<T,a,d> M((INITIAL_SIZE)Columns(),(INITIAL_SIZE)d);
        for(int i=0;i<Columns();i++)
            for(int k=0;k<Rows();k++)
                M(i,k)=(*this)(k,i)*A(k,k);
        return M;
    }

    template<class T_MATRIX0>
    auto Times_Transpose(const MATRIX_BASE<T,T_MATRIX0>& A) const
    {
        static const int a=MATRIX_INFO<T_MATRIX>::m;
        static const int b=MATRIX_INFO<T_MATRIX0>::m;
        SIZED_MATRIX<T,a,b> M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)A.Rows());
        Add_Times_Transpose(Derived(),A,M);
        return M;
    }

    template<int d>
    auto Times_Transpose(const DIAGONAL_MATRIX<T,d>& A) const
    {return Derived()*A;}

    template<int d>
    auto Times_Transpose(const SYMMETRIC_MATRIX<T,d>& A) const
    {return Derived()*A;}

    template<class T_VECTOR>
    auto Transpose_Times(const ARRAY_BASE<T,T_VECTOR>& y) const
    {
        assert(y.Size()==Rows());
        RIGHT_VECTOR result((INITIAL_SIZE)Columns());
        Add_Transpose_Times(Derived(),y.Derived(),result);
        return result;
    }

    template<class T_MATRIX0>
    auto operator*(const MATRIX_BASE<T,T_MATRIX0>& A) const
    {
        assert(Columns()==A.Rows());
        static const int a=MATRIX_INFO<T_MATRIX>::m;
        static const int b=MATRIX_INFO<T_MATRIX0>::n;
        SIZED_MATRIX<T,a,b> M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)A.Columns());
        Add_Times(Derived(),A,M);
        return M;
    }

    template<class T_VECTOR>
    auto operator*(const ARRAY_BASE<T,T_VECTOR>& y) const
    {
        assert(y.Size()==Columns());
        LEFT_VECTOR result((INITIAL_SIZE)Rows());
        Add_Times(Derived(),y.Derived(),result);
        return result;
    }

    template<class T_VECTOR,class T_VECTOR1>
    void Times(const ARRAY_BASE<T,T_VECTOR>& y,ARRAY_BASE<T,T_VECTOR1>& z) const
    {z.Fill(T());Add_Times(Derived(),y.Derived(),z.Derived());}

    template<int d>
    auto operator*(const DIAGONAL_MATRIX<T,d>& A) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);T_MATRIX M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int k=0;k<Columns();k++) for(int i=0;i<Rows();i++) M(i,k)=(*this)(i,k)*A(k,k);
    return M;}

    template<int d>
    auto operator*(const SYMMETRIC_MATRIX<T,d>& A) const
    {assert(Columns()==A.Rows());T_MATRIX M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)A.Columns());Add_Times(Derived(),A,M);return M;}

    template<int d>
    T_MATRIX operator*(const UPPER_TRIANGULAR_MATRIX<T,d>& A) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);assert(Columns()==d);
    T_MATRIX M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)d);
    for(int j=0;j<d;j++) for(int k=0;k<=j;k++) for(int i=0;i<Rows();i++) M(i,j)+=(*this)(i,k)*A(k,j);
    return M;}

    template<int d>
    T_MATRIX Times_Transpose(const UPPER_TRIANGULAR_MATRIX<T,d>& A) const
    {assert(Columns()==d);T_MATRIX M((INITIAL_SIZE)Rows(),(INITIAL_SIZE)d);
    for(int j=0;j<Rows();j++) for(int k=0;k<d;k++) for(int i=0;i<=k;i++) M(j,i)+=(*this)(j,k)*A(i,k);
    return M;}

    template<int d>
    auto Transpose_Times(const UPPER_TRIANGULAR_MATRIX<T,d>& A) const
    {
        assert(Rows()==d);
        static const int a=MATRIX_INFO<T_MATRIX>::n;
        SIZED_MATRIX<T,a,d> M((INITIAL_SIZE)Columns(),(INITIAL_SIZE)d);
        for(int j=0;j<Columns();j++)
            for(int k=0;k<d;k++)
                for(int i=0;i<=k;i++)
                    M(j,k)+=(*this)(i,j)*A(i,k);
        return M;
    }

    T_MATRIX& operator*=(const T a)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) (*this)(i,j)*=a;
    return Derived();}

    T_MATRIX& operator/=(const T a)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);
    return Derived()*=(1/a);}

    T_MATRIX operator*(const T a) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);T_MATRIX matrix((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) matrix(i,j)=(*this)(i,j)*a;
    return matrix;}

    T_MATRIX operator/(const T a) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);return Derived()*(1/a);}

    T_MATRIX& operator+=(const T a)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);assert(Rows()==Columns());
    for(int i=0;i<Rows();i++) (*this)(i,i)+=a;
    return Derived();}

    T_MATRIX& operator-=(const T a)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);assert(Rows()==Columns());
    for(int i=0;i<Rows();i++) (*this)(i,i)-=a;
    return Derived();}

    T_MATRIX operator+(const T a) const
    {return T_MATRIX(Derived())+=a;}

    T_MATRIX operator-(const T a) const
    {return T_MATRIX(Derived())-=a;}

    template<class T_MATRIX0>
    T_MATRIX& operator+=(const MATRIX_BASE<T,T_MATRIX0>& A)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,T_MATRIX0);assert(Rows()==A.Rows() && Columns()==A.Columns());
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) (*this)(i,j)+=A(i,j);
    return Derived();}

    template<int d>
    T_MATRIX& operator+=(const SYMMETRIC_MATRIX<T,d>& A)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,SYMMETRIC_MATRIX<T,d>);assert(Rows()==A.Rows() && Columns()==A.Columns());
    for(int j=0;j<Columns();j++){
        for(int i=j+1;i<Rows();i++){
            T element=A.Element_Lower(i,j);
            (*this)(i,j)+=element;
            (*this)(j,i)+=element;}
        (*this)(j,j)+=A.Element_Lower(j,j);}
    return Derived();}

    template<class T_MATRIX0>
    T_MATRIX& operator-=(const MATRIX_BASE<T,T_MATRIX0>& A)
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,T_MATRIX0);assert(Rows()==A.Rows() && Columns()==A.Columns());
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) (*this)(i,j)-=A(i,j);
    return Derived();}

    template<class T_MATRIX0>
    T_MATRIX operator+(const MATRIX_BASE<T,T_MATRIX0>& A) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,T_MATRIX0);assert(Rows()==A.Rows() && Columns()==A.Columns());T_MATRIX matrix((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) matrix(i,j)=(*this)(i,j)+A(i,j);
    return matrix;}

    template<class T_MATRIX0>
    T_MATRIX operator-(const MATRIX_BASE<T,T_MATRIX0>& A) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX,T_MATRIX0);assert(Rows()==A.Rows() && Columns()==A.Columns());T_MATRIX matrix((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) matrix(i,j)=(*this)(i,j)-A(i,j);
    return matrix;}

    T_MATRIX operator-() const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);T_MATRIX matrix((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) matrix(i,j)=-(*this)(i,j);
    return matrix;}

    T Trace() const
    {assert(Rows()==Columns());T trace=0;
    for(int i=0;i<Columns();i++) trace+=(*this)(i,i);
    return trace;}

    template<class T_VECTOR>
    T_MATRIX Permute_Columns(const ARRAY_BASE<int,T_VECTOR>& p) const
    {assert(Columns()==p.Size());T_MATRIX x((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int i=0;i<Rows();i++) for(int j=0;j<Columns();j++) x(i,j)=(*this)(i,p(j));
    return x;}

    template<class T_VECTOR>
    T_MATRIX Unpermute_Columns(const ARRAY_BASE<int,T_VECTOR>& p) const
    {assert(Columns()==p.Size());T_MATRIX x((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int i=0;i<Rows();i++) for(int j=0;j<Columns();j++) x(i,p(j))=(*this)(i,j);
    return x;}

    template<class T_VECTOR0,class T_VECTOR1>
    static T_MATRIX Outer_Product(const ARRAY_BASE<T,T_VECTOR0>& u,const ARRAY_BASE<T,T_VECTOR1>& v)
    {T_MATRIX result((INITIAL_SIZE)u.Size(),(INITIAL_SIZE)v.Size());
    for(int i=0;i<u.Size();i++) for(int j=0;j<v.Size();j++) result(i,j)=u(i)*v(j);
    return result;}

    template<class T_VECTOR0,class T_VECTOR1>
    void Gauss_Seidel_Single_Iteration(ARRAY_BASE<T,T_VECTOR0>& x,const ARRAY_BASE<T,T_VECTOR1>& b) const
    {assert(Rows()==Columns() && x.Size()==b.Size() && x.Size()==Columns());
    for(int i=0;i<Columns();i++){
        T rho=0;
        for(int j=0;j<i;j++) rho+=(*this)(i,j)*x(j);
        for(int j=i+1;j<Columns();j++) rho+=(*this)(i,j)*x(j);
        x(i)=(b(i)-rho)/(*this)(i,i);}}

    void Left_Givens_Rotation(const int i,const int j,const T c,const T s)
    {assert(0<=i && i<j && j<Rows());
    for(int k=0;k<Columns();k++){
        T x=(*this)(i,k);
        (*this)(i,k)=c*(*this)(i,k)-s*(*this)(j,k);
        (*this)(j,k)=s*x+c*(*this)(j,k);}}
    
    void Right_Givens_Rotation(const int i,const int j,const T c,const T s)
    {assert(0<=i && i<j && j<Columns());
    for(int k=0;k<Rows();k++){
        T x=(*this)(k,i);
        (*this)(k,i)=c*(*this)(k,i)-s*(*this)(k,j);
        (*this)(k,j)=s*x+c*(*this)(k,j);}}

    T Max_Abs() const
    {T max_abs=0;
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) max_abs=max(max_abs,abs((*this)(i,j)));
    return max_abs;}

    T Infinity_Norm() const
    {T max_sum=0;
    for(int j=0;j<Columns();j++){
        T sum=0;
        for(int i=0;i<Rows();i++) sum+=abs((*this)(i,j));
        max_sum=max(sum,max_sum);}
    return max_sum;}

    T Frobenius_Norm() const
    {return sqrt(Derived().Frobenius_Norm_Squared());}

    T Frobenius_Norm_Squared() const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);T sum=0;
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) sum+=sqr((*this)(i,j));
    return sum;}

    T Double_Contract(const MATRIX_BASE<T,T_MATRIX>& a) const
    {WARN_IF_NOT_EFFICIENT(T_MATRIX);T sum=0;
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) sum+=(*this)(i,j)*a(i,j);
    return sum;}

    template<int d>
    T Double_Contract(const SYMMETRIC_MATRIX<T,d>& a) const
    {assert(d==Rows() && d==Columns());T sum=0;
    for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) sum+=(*this)(i,j)*a(i,j);
    return sum;}

    template<class T_MATRIX1>
    void Set_Submatrix(const int istart,const int jstart,const MATRIX_BASE<T,T_MATRIX1>& a)
    {for(int i=0;i<a.Rows();i++) for(int j=0;j<a.Columns();j++) (*this)(istart+i,jstart+j)=a(i,j);}

    void Set_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,3>& a)
    {(*this)(istart,jstart)=a.x00;(*this)(istart,jstart+1)=a.x10;(*this)(istart+1,jstart)=a.x10;(*this)(istart,jstart+2)=a.x20;(*this)(istart+2,jstart)=a.x20;
    (*this)(istart+1,jstart+1)=a.x11;(*this)(istart+1,jstart+2)=a.x21;(*this)(istart+2,jstart+1)=a.x21;(*this)(istart+2,jstart+2)=a.x22;}

    void Set_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,2>& a)
    {(*this)(istart,jstart)=a.x00;(*this)(istart,jstart+1)=a.x10;(*this)(istart+1,jstart)=a.x10;(*this)(istart+1,jstart+1)=a.x11;}

    void Set_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,1>& a)
    {(*this)(istart,jstart)=a.x00;}

    void Set_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,0>& a)
    {}

    void Set_Submatrix(const int istart,const int jstart,const DIAGONAL_MATRIX<T,3>& a)
    {(*this)(istart,jstart)=a.x00;(*this)(istart+1,jstart+1)=a.x11;(*this)(istart+2,jstart+2)=a.x22;}

    void Set_Submatrix(const int istart,const int jstart,const DIAGONAL_MATRIX<T,2>& a)
    {(*this)(istart,jstart)=a.x00;(*this)(istart+1,jstart+1)=a.x11;}

    void Set_Submatrix(const int istart,const int jstart,const DIAGONAL_MATRIX<T,1>& a)
    {(*this)(istart,jstart)=a.x00;}

    void Set_Submatrix(const int istart,const int jstart,const DIAGONAL_MATRIX<T,0>& a)
    {}

    template<class T_VECTOR>
    void Set_Submatrix(const int istart,const int jstart,const ARRAY_BASE<T,T_VECTOR>& a)
    {for(int i=0;i<a.Size();i++) (*this)(istart+i,jstart)=a(i);}

    template<class T_MATRIX1>
    void Add_To_Submatrix(const int istart,const int jstart,const MATRIX_BASE<T,T_MATRIX1>& a)
    {for(int i=0;i<a.Rows();i++) for(int j=0;j<a.Columns();j++) (*this)(istart+i,jstart+j)+=a(i,j);}

    void Add_To_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,3>& a)
    {(*this)(istart,jstart)+=a.x00;(*this)(istart,jstart+1)+=a.x10;(*this)(istart+1,jstart)+=a.x10;(*this)(istart,jstart+2)+=a.x20;(*this)(istart+2,jstart)+=a.x20;
    (*this)(istart+1,jstart+1)+=a.x11;(*this)(istart+1,jstart+2)+=a.x21;(*this)(istart+2,jstart+1)+=a.x21;(*this)(istart+2,jstart+2)+=a.x22;}

    void Add_To_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,2>& a)
    {(*this)(istart,jstart)+=a.x00;(*this)(istart,jstart+1)+=a.x10;(*this)(istart+1,jstart)+=a.x10;(*this)(istart+1,jstart+1)+=a.x11;}

    void Add_To_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,1>& a)
    {(*this)(istart,jstart)+=a.x00;}

    void Add_To_Submatrix(const int istart,const int jstart,const SYMMETRIC_MATRIX<T,0>& a)
    {}

    template<int d>
    void Add_To_Submatrix(const int istart,const int jstart,const DIAGONAL_MATRIX<T,d>& a)
    {for(int i=0;i<d;i++) (*this)(istart+i,jstart+i)+=a.x(i);}

    template<class T_VECTOR>
    void Add_To_Submatrix(const int istart,const int jstart,const ARRAY_BASE<T,T_VECTOR>& a)
    {for(int i=0;i<a.Size();i++) (*this)(istart+i,jstart)+=a(i);}

    template<class T_MATRIX1>
    void Subtract_From_Submatrix(const int istart,const int jstart,const MATRIX_BASE<T,T_MATRIX1>& a)
    {for(int i=0;i<a.Rows();i++) for(int j=0;j<a.Columns();j++) (*this)(istart+i,jstart+j)-=a(i,j);}

    template<class T_MATRIX1>
    void Get_Submatrix(const int istart,const int jstart,MATRIX_BASE<T,T_MATRIX1>& a) const
    {for(int i=0;i<a.Rows();i++) for(int j=0;j<a.Columns();j++) a(i,j)=(*this)(istart+i,jstart+j);}

    void Get_Submatrix(const int istart,const int jstart,SYMMETRIC_MATRIX<T,3>& a) const
    {a.x00=(*this)(istart,jstart);a.x10=(*this)(istart,jstart+1);a.x10=(*this)(istart+1,jstart);a.x20=(*this)(istart,jstart+2);a.x20=(*this)(istart+2,jstart);
    a.x11=(*this)(istart+1,jstart+1);a.x21=(*this)(istart+1,jstart+2);a.x21=(*this)(istart+2,jstart+1);a.x22=(*this)(istart+2,jstart+2);}

    void Get_Submatrix(const int istart,const int jstart,DIAGONAL_MATRIX<T,3>& a) const
    {a.x00=(*this)(istart,jstart);a.x11=(*this)(istart+1,jstart+1);a.x22=(*this)(istart+2,jstart+2);}

    template<class T_VECTOR>
    void Set_Column(const int j,const ARRAY_BASE<T,T_VECTOR>& a)
    {for(int i=0;i<Rows();i++) (*this)(i,j)=a(i);}

    template<class T_VECTOR>
    void Add_Column(const int j,const ARRAY_BASE<T,T_VECTOR>& a)
    {for(int i=0;i<Rows();i++) (*this)(i,j)+=a(i);}

    template<class T_VECTOR>
    void Get_Column(const int j,ARRAY_BASE<T,T_VECTOR>& a) const
    {for(int i=0;i<Rows();i++) a(i)=(*this)(i,j);}

    template<class T_VECTOR>
    void Set_Row(const int i,const ARRAY_BASE<T,T_VECTOR>& a)
    {for(int j=0;j<Columns();j++) (*this)(i,j)=a(j);}

    template<class T_VECTOR>
    void Add_Row(const int i,const ARRAY_BASE<T,T_VECTOR>& a)
    {for(int j=0;j<Columns();j++) (*this)(i,j)+=a(j);}

    template<class T_VECTOR>
    void Get_Row(const int i,ARRAY_BASE<T,T_VECTOR>& a) const
    {for(int j=0;j<Columns();j++) a(j)=(*this)(i,j);}

    void Set_Zero_Matrix()
    {for(int i=0;i<Rows();i++) for(int j=0;j<Columns();j++) (*this)(i,j)=T();}

    void Add_Identity_Matrix()
    {assert(Rows()==Columns());
    for(int i=0;i<Columns();i++) (*this)(i,i)+=1;}

    void Set_Identity_Matrix()
    {Set_Zero_Matrix();Add_Identity_Matrix();}

    static T_MATRIX Identity_Matrix(const int n)
    {T_MATRIX A(n);A.Add_Identity_Matrix();return A;}

    template<class T_VECTOR>
    T Symmetric_Conjugate(const ARRAY_BASE<T,T_VECTOR>& v) const
    {assert(Rows()==Columns());T r=0;
    for(int j=0;j<Columns();j++){
        T a=0;
        for(int i=0;i<j;i++) a+=v(i)*(*this)(i,j);
        r+=(a+a+v(j)*(*this)(j,j))*v(j);}
    return r;}

    template<class T_VECTOR>
    RIGHT_VECTOR Lower_Triangular_Solve(const ARRAY_BASE<T,T_VECTOR>& b) const
    {assert(Rows()==Columns() && Columns()==b.Size());RIGHT_VECTOR x(INITIAL_SIZE(b.Size()));
    for(int i=0;i<Columns();i++){
        x(i)=b(i);
        for(int j=0;j<i;j++) x(i)-=(*this)(i,j)*x(j);
        x(i)/=(*this)(i,i);}
    return x;}

    template<class T_VECTOR>
    LEFT_VECTOR Transpose_Lower_Triangular_Solve(const ARRAY_BASE<T,T_VECTOR>& b) const
    {assert(Rows()==Columns() && Columns()==b.Size());LEFT_VECTOR x(b);
    for(int i=0;i<Columns();i++){
        x(i)/=(*this)(i,i);
        for(int j=i+1;j<Columns();j++) x(j)-=(*this)(i,j)*x(i);}
    return x;}

    template<class T_VECTOR>
    RIGHT_VECTOR Upper_Triangular_Solve(const ARRAY_BASE<T,T_VECTOR>& b) const
    {assert(Rows()==Columns() && Columns()==b.Size());RIGHT_VECTOR x(INITIAL_SIZE(b.Size()));
    for(int i=Columns()-1;i>=0;i--){
        x(i)=b(i);
        for(int j=Columns()-1;j>i;j--) x(i)-=(*this)(i,j)*x(j);
        x(i)/=(*this)(i,i);}
    return x;}

    template<class T_VECTOR>
    LEFT_VECTOR Transpose_Upper_Triangular_Solve(const ARRAY_BASE<T,T_VECTOR>& b) const
    {assert(Rows()==Columns() && Columns()==b.Size());LEFT_VECTOR x(b);
    for(int i=Columns()-1;i>=0;i--){
        x(i)/=(*this)(i,i);
        for(int j=i-1;j>=0;j--) x(j)-=(*this)(i,j)*x(i);}
    return x;}

    template<class T_MATRIX1>
    T_MATRIX1 Upper_Triangular_Solve(const MATRIX_BASE<T,T_MATRIX1>& b) const
    {assert(Rows()==Columns() && Columns()==b.Rows());T_MATRIX1 x(INITIAL_SIZE(b.Rows()),INITIAL_SIZE(b.Columns()));
    for(int bcol=0;bcol<b.Columns();bcol++)
        for(int i=Columns()-1;i>=0;i--){
            x(i,bcol)=b(i,bcol);
            for(int j=Columns()-1;j>=i+1;j--) x(i,bcol)-=(*this)(i,j)*x(j,bcol);
            x(i,bcol)/=(*this)(i,i);}
    return x;}

    template<class T_VECTOR>
    RIGHT_VECTOR In_Place_Cholesky_Solve(const ARRAY_BASE<T,T_VECTOR>& b)
    {assert(Rows()==Columns());In_Place_Cholesky_Factorization();return Transpose_Upper_Triangular_Solve(Lower_Triangular_Solve(b));}

    template<class T_VECTOR>
    RIGHT_VECTOR Cholesky_Solve(const ARRAY_BASE<T,T_VECTOR>& b) const
    {return T_MATRIX(Derived()).In_Place_Cholesky_Solve(b);}

    void In_Place_Cholesky_Inverse()
    {return T_MATRIX(Derived()).In_Place_Cholesky_Inverse(*this);}

    template<class T_MATRIX1>
    void Cholesky_Inverse(MATRIX_BASE<T,T_MATRIX1>& inverse) const
    {return T_MATRIX(Derived()).In_Place_Cholesky_Inverse(inverse);}

    template<class T_VECTOR>
    RIGHT_VECTOR In_Place_Gram_Schmidt_QR_Solve(const ARRAY_BASE<T,T_VECTOR>& b)
    {T_MATRIX R;In_Place_Gram_Schmidt_QR_Factorization(R);return R.Upper_Triangular_Solve(Derived().Transpose_Times(b.Derived()));}

    template<class T_VECTOR>
    RIGHT_VECTOR Gram_Schmidt_QR_Solve(const ARRAY_BASE<T,T_VECTOR>& b) const
    {return T_MATRIX(Derived()).In_Place_Gram_Schmidt_QR_Solve(b);}

    template<class T_VECTOR,class T_MATRIX1>
    static typename T_MATRIX1::LEFT_VECTOR Householder_Transform(const ARRAY_BASE<T,T_VECTOR>& b,const MATRIX_BASE<T,T_MATRIX1>& V)
    {assert(V.Rows()==b.Size());typename T_MATRIX1::LEFT_VECTOR result(b),v((INITIAL_SIZE)V.Rows());
    for(int j=0;j<V.Columns();j++){V.Get_Column(j,v);result=result.Householder_Transform(v);}
    return result;}

    template<class T_VECTOR>
    LEFT_VECTOR Householder_QR_Solve(const ARRAY_BASE<T,T_VECTOR>& b)
    {T_MATRIX V,R;Householder_QR_Factorization(V,R);LEFT_VECTOR c=Householder_Transform(b,V),c_short((INITIAL_SIZE)Columns());
    for(int i=0;i<Columns();i++) c_short(i)=c(i);
    return R.Upper_Triangular_Solve(c_short);}

    T Condition_Number() const
    {assert(Rows()==Columns());T_MATRIX inverse;PLU_Inverse(inverse);return Infinity_Norm()*inverse.Infinity_Norm();}

    template<class T_VECTOR>
    RIGHT_VECTOR In_Place_PLU_Solve(const ARRAY_BASE<T,T_VECTOR>& b)
    {assert(Rows()==Columns());COLUMN_PERMUTATION p;T_MATRIX L;In_Place_PLU_Factorization(L,p);
    return Upper_Triangular_Solve(L.Lower_Triangular_Solve(b.Subset(p)));}

    template<class T_VECTOR>
    RIGHT_VECTOR PLU_Solve(const ARRAY_BASE<T,T_VECTOR>& b) const
    {return T_MATRIX(Derived()).In_Place_PLU_Solve(b);}

    template<class T_MATRIX1>
    void PLU_Inverse(MATRIX_BASE<T,T_MATRIX1>& inverse) const
    {(T_MATRIX1(Derived())).In_Place_PLU_Inverse(inverse);}

    template<class T_VECTOR>
    RIGHT_VECTOR In_Place_LU_Solve(const ARRAY_BASE<T,T_VECTOR>& b)
    {assert(Rows()==Columns());T_MATRIX L;In_Place_LU_Factorization(L);return Upper_Triangular_Solve(L.Lower_Triangular_Solve(b));}

    template<class T_VECTOR>
    RIGHT_VECTOR LU_Solve(const ARRAY_BASE<T,T_VECTOR>& b) const
    {return T_MATRIX(Derived()).In_Place_LU_Inverse(b);}

    template<class T_MATRIX1>
    void LU_Inverse(MATRIX_BASE<T,T_MATRIX1>& inverse) const
    {(T_MATRIX1(Derived())).In_Place_LU_Inverse(inverse);}

    template<class T_VECTOR>
    RIGHT_VECTOR Inverse_Times(const ARRAY_BASE<T,T_VECTOR>& b)
    {return PLU_Solve(b);}

//#####################################################################
    template<class T_MATRIX1> void In_Place_Cholesky_Inverse(MATRIX_BASE<T,T_MATRIX1>& inverse);
    template<class T_MATRIX1> void In_Place_PLU_Inverse(MATRIX_BASE<T,T_MATRIX1>& inverse);
    template<class T_MATRIX1> void In_Place_LU_Inverse(MATRIX_BASE<T,T_MATRIX1>& inverse);
    template<class T_MATRIX1> void In_Place_Gram_Schmidt_QR_Factorization(MATRIX_BASE<T,T_MATRIX1>& R); // this=Q
    template<class T_MATRIX1,class T_MATRIX2> void Householder_QR_Factorization(MATRIX_BASE<T,T_MATRIX1>& V,MATRIX_BASE<T,T_MATRIX2>& R);
    template<class T_VECTOR0,class T_VECTOR1> void In_Place_Robust_Householder_QR_Solve(ARRAY_BASE<T,T_VECTOR0>& b,ARRAY_BASE<int,T_VECTOR1>& p); // this=Q
    template<class T_MATRIX1> void In_Place_PLU_Factorization(MATRIX_BASE<T,T_MATRIX1>& L,COLUMN_PERMUTATION& p); // this=U
    void In_Place_Cholesky_Factorization(); // this=L
    template<class T_MATRIX1> void In_Place_LU_Factorization(MATRIX_BASE<T,T_MATRIX1>& L); // this=U
    int Number_Of_Nonzero_Rows(const T threshold) const;
    void Jacobi_Singular_Value_Decomposition(ARRAY<VECTOR<int,2> >& left_givens_pairs,ARRAY<VECTOR<T,2> >& left_givens_coefficients,
        ARRAY<VECTOR<int,2> >& right_givens_pairs,ARRAY<VECTOR<T,2> >& right_givens_coefficients,const T tolerance=(T)1e-10,
        const int max_iterations=1000000);
//#####################################################################
};
template<class T,class T_MATRIX> inline std::ostream& operator<<(std::ostream& output_stream,const MATRIX_BASE<T,T_MATRIX>& A)
{
    output_stream<<"[";
    for(int i=0;i<A.Rows();i++){
        for(int j=0;j<A.Columns();j++){
            output_stream<<A(i,j);
            if(j<A.Columns()-1) output_stream<<" ";}
        if(i<A.Rows()-1) output_stream<<"; ";}
    output_stream<<"]";
    return output_stream;
}

template<class T,class T_MATRIX,int d>
auto operator*(const DIAGONAL_MATRIX<T,d>& A,const MATRIX_BASE<T,T_MATRIX>& B)
{
    assert(d==B.Rows());
    T_MATRIX M((INITIAL_SIZE)B.Rows(),(INITIAL_SIZE)B.Columns());
    for(int i=0;i<B.Rows();i++){
        T a=A(i,i);
        for(int k=0;k<B.Columns();k++) M(i,k)=a*B(i,k);}
    return M;
}

template<class T,class T_MATRIX,int d>
auto operator*(const SYMMETRIC_MATRIX<T,d>& A,const MATRIX_BASE<T,T_MATRIX>& B)
{T_MATRIX M((INITIAL_SIZE)B.Rows(),(INITIAL_SIZE)B.Columns());B.Add_Times(A,B.Derived(),M);return M;}

template<class T,class T_MATRIX>
T_MATRIX exp(const MATRIX_BASE<T,T_MATRIX>& M);
template<class T> class RANDOM_NUMBERS;
//#####################################################################
// Function Random_Fill_Uniform
//#####################################################################
template<class T,class T_MATRIX> void
Random_Fill_Uniform(RANDOM_NUMBERS<T>& rand,MATRIX_BASE<T,T_MATRIX>& m,const T a,const T b)
{
    for(int i=0;i<m.Rows();i++) for(int j=0;j<m.Columns();j++) m(i,j)=rand.Get_Uniform_Number(a,b);
}
}
#endif
