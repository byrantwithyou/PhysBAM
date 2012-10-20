#include "MATRIX_TESTS.h"
using namespace PhysBAM;
//#####################################################################
// Function Arbitrary_Test_One_Size
//#####################################################################
template<class T> template<class T_MATRIX> bool MATRIX_TESTS<T>::
Arbitrary_Test_One_Size(T_MATRIX A,int count) const
{
    T_MATRIX B((INITIAL_SIZE)A.Rows(),(INITIAL_SIZE)A.Columns()),C(A);bool ok=true;T tolerance=std::numeric_limits<T>::epsilon();
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
//#####################################################################
// Function Arbitrary_Test_Two_Sizes
//#####################################################################
template<class T> template<class T_MATRIX1,class T_MATRIX2> bool MATRIX_TESTS<T>::
Arbitrary_Test_Two_Sizes(T_MATRIX1 A,T_MATRIX2 B,int count) const
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
//#####################################################################
// Function Test
//#####################################################################
template<class T> bool MATRIX_TESTS<T>::
Test(const bool test,const std::string& message,bool& ok) const
{
    if(test) return true;
    LOG::cout<<"FAILED: "<<message<<std::endl;
    ok=false;
    return false;
}
//#####################################################################
// Function Test
//#####################################################################
template<class T> template<class T_MATRIX> bool MATRIX_TESTS<T>::
Test(const bool test,const std::string& message,bool& ok,const T_MATRIX& A) const
{
    if(Test(test,message,ok)) return true;
    Identify_Matrix(A);
    return false;
}
//#####################################################################
// Function Test
//#####################################################################
template<class T> template<class T_MATRIX1,class T_MATRIX2> bool MATRIX_TESTS<T>::
Test(const bool test,const std::string& message,bool& ok,const T_MATRIX1& A,const T_MATRIX2& B) const
{
    if(Test(test,message,ok,A)) return true;
    Identify_Matrix(B);
    return false;
}
//#####################################################################
// Function Test
//#####################################################################
template<class T> template<class T_MATRIX1,class T_MATRIX2,class T_MATRIX3> bool MATRIX_TESTS<T>::
Test(const bool test,const std::string& message,bool& ok,const T_MATRIX1& A,const T_MATRIX2& B,const T_MATRIX3& C) const
{
    if(Test(test,message,ok,A,B)) return true;
    Identify_Matrix(C);
    return false;
}
//#####################################################################
// Function Conversion_Test
//#####################################################################
template<class T> template<class T_MATRIX> void MATRIX_TESTS<T>::
Conversion_Test(T_MATRIX& A,MATRIX_MXN<T>& B,bool& ok) const
{
    Test(Assert_Equal(MATRIX_MXN<T>(T_MATRIX(B)),B),"Conversion both ways is exact.",ok,A);
}
//#####################################################################
// Function Conversion_Test
//#####################################################################
template<class T> template<int d> void MATRIX_TESTS<T>::
Conversion_Test(DIAGONAL_MATRIX<T,d>& A,MATRIX_MXN<T>& B,bool& ok) const
{
}
//#####################################################################
// Function Conversion_Test
//#####################################################################
template<class T> template<int d> void MATRIX_TESTS<T>::
Conversion_Test(SYMMETRIC_MATRIX<T,d>& A,MATRIX_MXN<T>& B,bool& ok) const
{
}
//#####################################################################
// Function Identify_Matrix
//#####################################################################
template<class T> template<class T_MATRIX> void MATRIX_TESTS<T>::
Identify_Matrix(const T_MATRIX& a) const
{
    LOG::cout<<"Size: "<<a.Rows()<<" x "<<a.Columns()<<"   Type: "<<typeid(a).name()<<std::endl;
}
