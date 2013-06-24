//#####################################################################
// Copyright 2007-2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MATRIX_TESTS__
#define __MATRIX_TESTS__
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_0X0.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_1X1.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Utilities/TEST_BASE.h>
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
template<class T>
class MATRIX_TESTS:public TEST_BASE
{
    typedef VECTOR<T,3> TV;
    mutable RANDOM_NUMBERS<T> rand;
public:

    MATRIX_TESTS()
        :TEST_BASE("matrix")
    {}

    virtual ~MATRIX_TESTS(){}

    TEST_RESULT Run_Test(int n) PHYSBAM_OVERRIDE;
    int Number_Of_Tests() const PHYSBAM_OVERRIDE;
    template<class T_MATRIX> bool Assert_Zero(const T_MATRIX& a,T tolerance=0) const;
    template<class T_MATRIX1,class T_MATRIX2> bool Assert_Equal(const T_MATRIX1& a,const T_MATRIX2& b,T tolerance=0) const;
    template<class T_MATRIX> void Identify_Matrix(const T_MATRIX& a) const;
    bool Test(const bool test,const std::string& message,bool& ok) const;
    template<class T_MATRIX1> bool Test(const bool test,const std::string& message,bool& ok,const T_MATRIX1& A) const;
    template<class T_MATRIX1,class T_MATRIX2> bool Test(const bool test,const std::string& message,bool& ok,const T_MATRIX1& A,const T_MATRIX2& B) const;
    template<class T_MATRIX1,class T_MATRIX2,class T_MATRIX3> bool Test(const bool test,const std::string& message,bool& ok,const T_MATRIX1& A,const T_MATRIX2& B,const T_MATRIX3& C) const;
    bool Dynamic_Tests_One_Size(int m,int n) const;
    bool Dynamic_Tests_Two_Sizes(int m,int n,int p) const;
    bool Dynamic_Tests_Three_Sizes(int m,int n,int p,int q) const;
    bool Dynamic_Tests(int size,int count) const;
    template<class T_MATRIX> void Conversion_Test(T_MATRIX& A,MATRIX_MXN<T>& B,bool& ok) const;
    template<int d> void Conversion_Test(DIAGONAL_MATRIX<T,d>& A,MATRIX_MXN<T>& B,bool& ok) const;
    template<int d> void Conversion_Test(SYMMETRIC_MATRIX<T,d>& A,MATRIX_MXN<T>& B,bool& ok) const;
    template<int d> void Conversion_Test(UPPER_TRIANGULAR_MATRIX<T,d>& A,MATRIX_MXN<T>& B,bool& ok) const;
    template<class T_MATRIX1> bool Arbitrary_Test_One_Size(T_MATRIX1 A,int count) const;
    template<class T_MATRIX1,class T_MATRIX2> bool Arbitrary_Test_Two_Sizes(T_MATRIX1 A,T_MATRIX2 B,int count) const;
    template<class T_MATRIX1,class T_MATRIX2> bool Arbitrary_Test_Two_Sizes_XT(T_MATRIX1 A,T_MATRIX2 B,int count) const;
    template<class T_MATRIX1,class T_MATRIX2> bool Arbitrary_Test_Two_Sizes_TX(T_MATRIX1 A,T_MATRIX2 B,int count) const;
//#####################################################################
};
}
#endif
