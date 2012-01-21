//#####################################################################
// Copyright 2007-2008, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_TESTS
//#####################################################################
//   1. Arithematic tests.
//#####################################################################

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/TEST_BASE.h>
#include <limits>
namespace PhysBAM{

/*
NOT YET TESTED:
Anything with UPPER_TRIANGULAR_MATRIX
Anything with Cholesky LU PLU QR
Anything with vectors
Anything with Cross_Product_Matrix
Set_Submatrix Get_Submatrix Add_To_Submatrix Subtract_From_Submatrix
Set_Column Get_Column Set_Row Get_Row
Max_Abs Infinity_Norm Frobenius_Norm Frobenius_Norm_Squared Parallelepiped_Measure Trace
Upper_Triangular_Solve Condition_Number Number_Of_Nonzero_Rows
Column Permute_Columns Unpermute_Columns Outer_Product
Normal_Equations_Matrix Normal_Equations_Solve Solve_Linear_System
Column_Sum Componentwise_Min Componentwise_Max Determinant
Higham_Iterate Fast_Singular_Value_Decomposition Fast_Indefinite_Polar_Decomposition

TESTED ONLY FOR DYNAMIC SIZES:
matrix+=matrix matrix-=matrix
Set_Zero_Matrix Add_Identity_Matrix Set_Identity_Matrix
matrix==matrix matrix!=matrix matrix+scalar matrix+=scalar
*/

template<class T_input>
class MATRIX_TESTS:public TEST_BASE
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;
    mutable RANDOM_NUMBERS<T> rand;
public:

    MATRIX_TESTS()
        :TEST_BASE("matrix")
    {}

    virtual ~MATRIX_TESTS(){}

    TEST_RESULT Run_Test(int n) PHYSBAM_OVERRIDE
    {switch(n){
            case 1:if(Dynamic_Tests(6,20)) return success;return failure;
            case 2:if(TEST_GENERATOR<T,1,6,6*6*6-1>(*this,20).v) return success;return failure;}
    return failure;}

    int Number_Of_Tests() const PHYSBAM_OVERRIDE
    {return 2;}

    template<class T_MATRIX>
    bool Assert_Zero(const T_MATRIX& a,T tolerance=0) const
    {if(a.Max_Abs()>tolerance) LOG::cout<<a.Max_Abs()<<"  "<<tolerance<<std::endl;return a.Max_Abs()<=tolerance;}

    template<class T_MATRIX1,class T_MATRIX2>
    bool Assert_Equal(const T_MATRIX1& a,const T_MATRIX2& b,T tolerance=0) const
    {tolerance*=max(a.Max_Abs(),b.Max_Abs(),(T)1);if(!Assert_Zero(a-b,tolerance)) LOG::cout<<a<<std::endl<<b<<std::endl;return Assert_Zero(a-b,tolerance);}

    template<class T_MATRIX>
    void Identify_Matrix(const T_MATRIX& a) const
    {LOG::cout<<"Size: "<<a.Rows()<<" x "<<a.Columns()<<"   Type: "<<typeid(a).name()<<std::endl;}

    bool Test(const bool test,const std::string& message,bool& ok) const
    {if(test) return true;LOG::cout<<"FAILED: "<<message<<std::endl;ok=false;return false;}

    template<class T_MATRIX1>
    bool Test(const bool test,const std::string& message,bool& ok,const T_MATRIX1& A) const
    {if(Test(test,message,ok)) return true;Identify_Matrix(A);return false;}

    template<class T_MATRIX1,class T_MATRIX2>
    bool Test(const bool test,const std::string& message,bool& ok,const T_MATRIX1& A,const T_MATRIX2& B) const
    {if(Test(test,message,ok,A)) return true;Identify_Matrix(B);return false;}

    template<class T_MATRIX1,class T_MATRIX2,class T_MATRIX3>
    bool Test(const bool test,const std::string& message,bool& ok,const T_MATRIX1& A,const T_MATRIX2& B,const T_MATRIX3& C) const
    {if(Test(test,message,ok,A,B)) return true;Identify_Matrix(C);return false;}

    bool Dynamic_Tests_One_Size(int m,int n) const
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

    bool Dynamic_Tests_Two_Sizes(int m,int n,int p) const
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

    bool Dynamic_Tests_Three_Sizes(int m,int n,int p,int q) const
    {
        MATRIX_MXN<T> A(m,n),B(n,p),C(p,q);rand.Fill_Uniform(A,-1,1);rand.Fill_Uniform(B,-1,1);rand.Fill_Uniform(C,-1,1);
        T tolerance=std::numeric_limits<T>::epsilon();bool ok=true;

        Test(Assert_Equal((A*B)*C,A*(B*C),tolerance*B.Rows()*B.Columns()*2),"Multiplication associates.",ok,A,B,C);
        return ok;
    }

    bool Dynamic_Tests(int size,int count) const
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

    template<class T,class T_MATRIX>
    void Conversion_Test(T_MATRIX& A,MATRIX_MXN<T>& B,bool& ok) const
    {
        Test(Assert_Equal(MATRIX_MXN<T>(T_MATRIX(B)),B),"Conversion both ways is exact.",ok,A);
    }

    template<class T,int d>
    void Conversion_Test(DIAGONAL_MATRIX<T,d>& A,MATRIX_MXN<T>& B,bool& ok) const
    {}

    template<class T,int d>
    void Conversion_Test(SYMMETRIC_MATRIX<T,d>& A,MATRIX_MXN<T>& B,bool& ok) const
    {}

    template<class T_MATRIX1>
    bool Arbitrary_Test_One_Size(T_MATRIX1 A,int count) const
    {
        T_MATRIX1 B((INITIAL_SIZE)A.Rows(),(INITIAL_SIZE)A.Columns()),C(A);bool ok=true;T tolerance=std::numeric_limits<T>::epsilon();
        Test(B.Columns()==A.Columns() && B.Rows()==A.Rows(),"Dimension tests (B).",ok,A);
        Test(C.Columns()==A.Columns() && C.Rows()==A.Rows(),"Dimension tests (C).",ok,A);

        for(int i=0;i<count;i++){
            MATRIX_MXN<T> D(A.Rows(),A.Columns()),E(A.Rows(),A.Columns()),F;rand.Fill_Uniform(A,-1,1);rand.Fill_Uniform(B,-1,1);D=A;E=B;
            T s=rand.Get_Uniform_Number(-1,1);

            Conversion_Test(A,D,ok);

            Test(Assert_Equal(MATRIX_MXN<T>(-A),-D,tolerance),"Negation matches.",ok,A);
            Test(Assert_Equal(MATRIX_MXN<T>(s*A),s*D,tolerance),"Left scaling matches.",ok,A);
            Test(Assert_Equal(MATRIX_MXN<T>(A*s),D*s,tolerance),"Right scaling matches.",ok,A);
            Test(Assert_Equal(MATRIX_MXN<T>(A/s),D/s,tolerance),"Scalar division matches.",ok,A);
            Test(Assert_Equal(MATRIX_MXN<T>(A+B),D+E,tolerance),"Addition matches.",ok,A);
            Test(Assert_Equal(MATRIX_MXN<T>(A-B),D-E,tolerance),"Subtraction matches.",ok,A);
            Test(Assert_Equal(MATRIX_MXN<T>(A.Transposed()),D.Transposed(),tolerance),"Transpsoed matches.",ok,A);

            C=A;Test(Assert_Equal(MATRIX_MXN<T>(C),D,tolerance),"Assignment matches.",ok,A);
            C=A;C*=s;Test(Assert_Equal(MATRIX_MXN<T>(C),D*s,tolerance),"Scalar times equal matches.",ok,A);
            C=A;C/=s;Test(Assert_Equal(MATRIX_MXN<T>(C),D/s,tolerance),"Scalar divide equal matches.",ok,A);
            C=A;C+=B;Test(Assert_Equal(MATRIX_MXN<T>(C),D+E,tolerance),"Plus equal matches.",ok,A);
            C=A;C-=B;Test(Assert_Equal(MATRIX_MXN<T>(C),D-E,tolerance),"Minus equal matches.",ok,A);

            Test(Assert_Equal(MATRIX_MXN<T>(A),D),"Inputs not changed.",ok,A);
            Test(Assert_Equal(MATRIX_MXN<T>(B),E),"Inputs not changed.",ok,A);}
        return ok;
    }

    template<class T_MATRIX1,class T_MATRIX2>
    bool Arbitrary_Test_Two_Sizes(T_MATRIX1 A,T_MATRIX2 B,int count) const
    {
        T tolerance=std::numeric_limits<T>::epsilon()*2;bool ok=true;
        Test(A.Columns()==B.Rows(),"Dimension tests (A).",ok,A);

        for(int i=0;i<count && ok;i++){
            MATRIX_MXN<T> D(A.Rows(),A.Columns()),E(B.Rows(),B.Columns());rand.Fill_Uniform(A,-1,1);rand.Fill_Uniform(B,-1,1);D=MATRIX_MXN<T>(A);E=MATRIX_MXN<T>(B);

            Test(Assert_Equal(MATRIX_MXN<T>(A*B),D*E,tolerance),"Multiplication matches.",ok,A,B);
            Test(Assert_Equal(MATRIX_MXN<T>(A.Transposed().Transpose_Times(B)),D*E,tolerance),"Transpose_Times matches.",ok,A,B);
            Test(Assert_Equal(MATRIX_MXN<T>(A.Times_Transpose(B.Transposed())),D*E,tolerance),"Times_Transpose matches.",ok,A,B);

            Test(Assert_Equal(MATRIX_MXN<T>(A),D),"Inputs not changed.",ok,A,B);
            Test(Assert_Equal(MATRIX_MXN<T>(B),E),"Inputs not changed.",ok,A,B);}
        return ok;
    }

    template<class T2,class M1,class M2,class M3,int a,int b,int c>
    struct TEST_GENERATOR2
    {bool v;TEST_GENERATOR2(const MATRIX_TESTS<T2>& m,int n){v=regular(m,n)&square(m,n,VECTOR<int,a==b && (a==2 || a==3)>());}
    bool regular(const MATRIX_TESTS<T2>& m,int n)
        {TEST_GENERATOR2<T2,M2,M3,MATRIX<T2,a,b>,b,c,a> tg1(m,n);TEST_GENERATOR2<T2,M2,M3,MATRIX_MXN<T2>,b,c,a> tg2(m,n);return tg1.v && tg2.v;}
    bool square(const MATRIX_TESTS<T2>& m,int n,VECTOR<int,0> w){return true;}
    bool square(const MATRIX_TESTS<T2>& m,int n,VECTOR<int,1> w)
        {TEST_GENERATOR2<T2,M2,M3,SYMMETRIC_MATRIX<T2,a>,b,c,a> tg1(m,n);TEST_GENERATOR2<T2,M2,M3,DIAGONAL_MATRIX<T2,a>,b,c,a> tg2(m,n);return tg1.v && tg2.v;}};

    template<class T2,class M2,class M3,int a,int b,int c>
    struct TEST_GENERATOR2<T2,bool,M2,M3,a,b,c>
    {bool v;TEST_GENERATOR2(const MATRIX_TESTS<T2>& m,int n){v=true;
        for(int i=0;i<n;i++) v=m.Arbitrary_Test_Two_Sizes(M2(INITIAL_SIZE(b),INITIAL_SIZE(c)),M3(INITIAL_SIZE(c),INITIAL_SIZE(a)),n)&v;
        if(a==b) for(int i=0;i<n;i++) v=m.Arbitrary_Test_One_Size(M2(INITIAL_SIZE(b),INITIAL_SIZE(c)),n)&v;}};

    template<class T2,int a,int b,int x> // use x=b*b*b*b-1
    struct TEST_GENERATOR
    {bool v;TEST_GENERATOR(const MATRIX_TESTS<T2>& m,int n){v=TEST_GENERATOR2<T2,int,int,bool,x/b/b+a,x/b%b+a,x%b+a>(m,n).v;v=(TEST_GENERATOR<T2,a,b,x-1>(m,n).v && v);}};

    template<class T2,int a,int b>
    struct TEST_GENERATOR<T2,a,b,-1>
    {bool v;TEST_GENERATOR(const MATRIX_TESTS<T2>& m,int n){v=true;}};
//#####################################################################
};
}
using namespace PhysBAM;
#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
MATRIX_TESTS<float> matrix_tests;
#else
MATRIX_TESTS<double> matrix_tests;
#endif
