//#####################################################################
// Copyright 2007-2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TENSOR_TESTS__
#define __TENSOR_TESTS__
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

template<class T>
class TENSOR_TESTS:public TEST_BASE
{
    typedef VECTOR<T,3> TV;
public:
    mutable RANDOM_NUMBERS<T> rand;
    mutable bool status;

    TENSOR_TESTS()
        :TEST_BASE("tensor")
    {}

    virtual ~TENSOR_TESTS(){}

    TEST_RESULT Run_Test(int n) override;
    int Number_Of_Tests() const override;
    template<class T_VECTOR> bool Assert_Zero_Vector(const T_VECTOR& a,const char* str,T tolerance) const;
    template<class T_VECTOR0,class T_VECTOR1> bool Assert_Equal_Vector(const T_VECTOR0& a,const T_VECTOR1& b,const char* str,T tolerance) const;
    template<class T_MATRIX> bool Assert_Zero_Matrix(const T_MATRIX& a,const char* str,T tolerance) const;
    template<class T_MATRIX0,class T_MATRIX1> bool Assert_Equal_Matrix(const T_MATRIX0& a,const T_MATRIX1& b,const char* str,T tolerance) const;
    template<class T_TENSOR> bool Assert_Zero_Tensor(const T_TENSOR& a,const char* str,T tolerance) const;
    template<class T_TENSOR0,class T_TENSOR1> bool Assert_Equal_Tensor(const T_TENSOR0& a,const T_TENSOR1& b,const char* str,T tolerance) const;

    template<int m,int n,int p> TENSOR<T,m,n,p> Fill_gt();
    template<int m,int n,int p> SYMMETRIC_TENSOR<T,m,n,p> Fill_st();
    template<int m,int n,int p> ZERO_TENSOR<T,m,n,p> Fill_zt();
    template<int u,int m,int n> VEC_ID_TENSOR<T,u,m,n> Fill_vt();
    template<int u,int m> VEC_ID_SYM_TENSOR<T,u,m> Fill_it();
    PERMUTATION_TENSOR<T> Fill_pt();

    template<int m,int n> MATRIX<T,m,n> Fill_gm();
    template<int m> SYMMETRIC_MATRIX<T,m> Fill_sm();
    template<int m,int n> ZERO_MATRIX<T,m,n> Fill_zm();
    template<int m> SCALE_MATRIX<T,m> Fill_cm();
    template<int m> IDENTITY_MATRIX<T,m> Fill_im();

    template<int m> VECTOR<T,m> Fill_gv();
    template<int m> ZERO_VECTOR<T,m> Fill_zv();

    void Generate_Tests();
    template<class T_TENSOR> void Test_Contract(const T_TENSOR& t);
    template<int s,class T_TENSOR,class T_VECTOR> void Test_Contract_TV(const T_TENSOR& t,const T_VECTOR& v);
    template<int r,int s,class T_TENSOR> void Test_Contract_T(const T_TENSOR& t);
    template<int r,int s,class T_TENSOR,class T_MATRIX> void Test_Contract_TM(const T_TENSOR& t,const T_MATRIX& m);
    template<int s,class T_TENSOR,class T_VECTOR> void Test_Tensor_Product(const T_TENSOR& m,const T_VECTOR& v);
    template<int r,int s,class T_TENSOR> void Test_Transposed(const T_TENSOR& m);
    template<int r,int s,class T_TENSOR> void Test_Twice_Symmetric_Part(const T_TENSOR& m);
    template<class T_TENSOR0,class T_TENSOR1> void Test_Add_Sub(const T_TENSOR0& a,const T_TENSOR1& b);

    void Test_0();
    template<int m> void Test_1();
    template<int m,int n> void Test_2();
    template<int m,int n,int p> void Test_3();
    template<int m,int n,int p,int q> void Test_4();

//#####################################################################
};
}
#endif
