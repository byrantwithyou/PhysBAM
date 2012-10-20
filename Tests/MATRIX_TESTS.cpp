//#####################################################################
// Copyright 2007-2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "MATRIX_TESTS.h"
#include "MATRIX_TEST_GENERATOR.h"
using namespace PhysBAM;
//#####################################################################
// Function Run_Test
//#####################################################################
template<class T> typename MATRIX_TESTS<T>::TEST_RESULT MATRIX_TESTS<T>::
Run_Test(int n) PHYSBAM_OVERRIDE
{
    switch(n){
        case 0:if(Dynamic_Tests(6,20)) return success;return failure;
        case 1:if(TEST_GENERATOR<T,1,6,6*6*6-1>(*this,20).v) return success;return failure;}
    return failure;
}
//#####################################################################
// Function Number_Of_Tests
//#####################################################################
template<class T> int MATRIX_TESTS<T>::
Number_Of_Tests() const PHYSBAM_OVERRIDE
{
    return 2;
}
//#####################################################################
// Function Assert_Zero
//#####################################################################
template<class T> template<class T_MATRIX> bool MATRIX_TESTS<T>::
Assert_Zero(const T_MATRIX& a,T tolerance) const
{
    if(a.Max_Abs()>tolerance) LOG::cout<<a.Max_Abs()<<"  "<<tolerance<<std::endl;
    return a.Max_Abs()<=tolerance;
}
//#####################################################################
// Function Assert_Equal
//#####################################################################
template<class T> template<class T_MATRIX1,class T_MATRIX2> bool MATRIX_TESTS<T>::
Assert_Equal(const T_MATRIX1& a,const T_MATRIX2& b,T tolerance) const
{
    tolerance*=max(a.Max_Abs(),b.Max_Abs(),(T)1);
    if(!Assert_Zero(a-b,tolerance)) LOG::cout<<a<<std::endl<<b<<std::endl;
    return Assert_Zero(a-b,tolerance);
}
//#####################################################################
// Function Dynamic_Tests_One_Size
//#####################################################################
template<class T> bool MATRIX_TESTS<T>::
Dynamic_Tests_One_Size(int m,int n) const
{
    MATRIX_MXN<T> A(m,n),B(A),C,D(A);rand.Fill_Uniform(A,-1,1);rand.Fill_Uniform(B,-1,1);rand.Fill_Uniform(D,-1,1);
    T s=rand.Get_Uniform_Number(-1,1),t=rand.Get_Uniform_Number(-1,1);bool ok=true;
    T tolerance=std::numeric_limits<T>::epsilon();

    Test(A.Rows()==m && A.Columns()==n,"Dimension tests.",ok,A);
    Test(B.Rows()==m && B.Columns()==n,"Dimension tests.",ok,A);

    Test(Assert_Zero(A-A),"Subtraction with self is zero.",ok,A);
    Test(Assert_Equal(- -A,A),"Negation is its own inverse.",ok,A);
    Test(Assert_Equal(A+A+A,A*3,tolerance*2),"Integer scaling as addition.",ok,A);
    Test(Assert_Equal(-A,A*-1),"Integer scaling as negation.",ok,A);
    Test(Assert_Equal(A*s,s*A,tolerance),"Scalar multiplication commutes.",ok,A);
    Test(Assert_Equal(A*(1/s),A*(1/s),tolerance),"Scalar division.",ok,A);
    Test(Assert_Equal(s*(t*A),(s*t)*A,tolerance*2),"Scalar multiplication associates.",ok,A);

    Test(Assert_Equal(A+B,B+A),"Addition commutes.",ok,A);
    Test(Assert_Equal(A-B,A+-B),"Definition of subtraction.",ok,A);
    Test(Assert_Equal(-(A-B),B-A),"Subtraction anticommutes.",ok,A);
    Test(Assert_Equal((B+A)+D,B+(A+D),tolerance*2),"Addition associates.",ok,A);
    Test(Assert_Equal(s*A+s*B,s*(A+B),tolerance*2),"Distributivity of scalar multiplication and addition.",ok,A);

    Test(Assert_Equal(A.Transposed().Transposed(),A),"Transpose is its own inverse.",ok,A);
    Test(Assert_Equal(A.Transposed()+B.Transposed(),(A+B).Transposed()),"Transpose of sum.",ok,A);
    C=A;C.Transpose();Test(Assert_Equal(C,A.Transposed()),"Transpose vs Transposed.",ok,A);

    C=A;Test(Assert_Equal(A,C),"Assignment.",ok,A);
    C=A;C+=B;Test(Assert_Equal(A+B,C,tolerance),"Plus equals.",ok,A);
    C=A;C-=B;Test(Assert_Equal(A-B,C,tolerance),"Minus equals.",ok,A);
    C=A;C*=s;Test(Assert_Equal(A*s,C,tolerance),"Scalar times equals.",ok,A);
    C=A;C/=s;Test(Assert_Equal(A/s,C,tolerance),"Scalar divide equals.",ok,A);

    Test(A==A,"Equality on equal matrices.",ok,A);
    if(m*n) Test(!(A==B),"Equality on different matrices.",ok,A);
    Test(!(A!=A),"Inequality on equal matrices.",ok,A);
    if(m*n) Test(A!=B,"Inequality on different matrices.",ok,A);

    Test(Assert_Equal(A.Identity_Matrix(m)*B,B,tolerance),"Left multiply identity.",ok,A);
    Test(Assert_Equal(A*A.Identity_Matrix(n),A,tolerance),"Right multiply identity.",ok,A);
    if(m==n){
        C=A;C.Add_Identity_Matrix();Test(Assert_Equal(C,A+A.Identity_Matrix(m),tolerance),"Add_Identity_Matrix.",ok,A);
        C.Set_Identity_Matrix();Test(Assert_Equal(A.Identity_Matrix(m),C,tolerance),"Set_Identity_Matrix.",ok,A);
        Test(Assert_Equal(A+s,A+s*A.Identity_Matrix(m),tolerance),"Addition with scalar.",ok,A);
        C=A;C+=s;Test(Assert_Equal(C,A+s*A.Identity_Matrix(m),tolerance),"Plus equals with scalar.",ok,A);}

    return ok;
}
//#####################################################################
// Function Dynamic_Tests_Two_Sizes
//#####################################################################
template<class T> bool MATRIX_TESTS<T>::
Dynamic_Tests_Two_Sizes(int m,int n,int p) const
{
    MATRIX_MXN<T> A(m,n),B(m,n),C(n,p),D(n,p);rand.Fill_Uniform(A,-1,1);rand.Fill_Uniform(B,-1,1);rand.Fill_Uniform(C,-1,1);
    rand.Fill_Uniform(D,-1,1);T s=rand.Get_Uniform_Number(-1,1),tolerance=std::numeric_limits<T>::epsilon();bool ok=true;

    Test(Assert_Equal((A+B)*C,A*C+B*C,tolerance*A.Columns()*2),"Right distributivity.",ok,A,C);
    Test(Assert_Equal(A*(C+D),A*C+A*D,tolerance*A.Columns()*2),"Left distributivity.",ok,A,C);
    Test(Assert_Equal((s*A)*C,s*(A*C),tolerance*A.Columns()*2),"Scalar and matrix multiplication associate.",ok,A,C);

    Test(Assert_Equal(C.Transposed()*A.Transposed(),(A*C).Transposed()),"Transpose of product.",ok,A,C);
    Test(Assert_Equal(A*C,A.Transposed().Transpose_Times(C)),"Transpose_Times.",ok,A,C);
    Test(Assert_Equal(A*C,A.Times_Transpose(C.Transposed())),"Times_Transpose.",ok,A,C);
    return ok;
}
//#####################################################################
// Function Dynamic_Tests_Three_Sizes
//#####################################################################
template<class T> bool MATRIX_TESTS<T>::
Dynamic_Tests_Three_Sizes(int m,int n,int p,int q) const
{
    MATRIX_MXN<T> A(m,n),B(n,p),C(p,q);rand.Fill_Uniform(A,-1,1);rand.Fill_Uniform(B,-1,1);rand.Fill_Uniform(C,-1,1);
    T tolerance=std::numeric_limits<T>::epsilon();bool ok=true;

    Test(Assert_Equal((A*B)*C,A*(B*C),tolerance*B.Rows()*B.Columns()*2),"Multiplication associates.",ok,A,B,C);
    return ok;
}
//#####################################################################
// Function Dynamic_Tests
//#####################################################################
template<class T> bool MATRIX_TESTS<T>::
Dynamic_Tests(int size,int count) const
{
    bool ok=true;
    for(int c=0;c<count;c++)
        for(int i=0;i<=size;i++) for(int j=0;j<size;j++){
                ok=Dynamic_Tests_One_Size(i,j)&ok;
                for(int k=0;k<=size;k++){
                    ok=Dynamic_Tests_Two_Sizes(i,j,k)&ok;
                    for(int m=0;m<=size;m++) ok=Dynamic_Tests_Three_Sizes(i,j,k,m)&ok;}}

    return ok;
}
#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
MATRIX_TESTS<float> matrix_tests;
#else
MATRIX_TESTS<double> matrix_tests;
#endif
