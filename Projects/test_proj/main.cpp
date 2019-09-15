//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Matrices/LAPACK.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>

using namespace PhysBAM;

template<class T,class A,class B,class C>
void Test_V_4(A& M,const B& b,C& x)
{
    LOG::printf("try: %s %s %s\n",typeid(A).name(),typeid(B).name(),typeid(C).name());
    MATRIX_MXN<T> M2(M);
    ARRAY<T> x2(x);
    Least_Squares_Solve(M,b,x);
    LOG::printf("sol: %P   %P\n",(M2.Transpose_Times(M2*x-b)).Magnitude(),x);
    M=M2;
    x=x2;
    Least_Squares_Solve_Raw(M,b,x);
    LOG::printf("row: %P   %P\n",(M2.Transpose_Times(M2*x-b)).Magnitude(),x);
    M=M2;
    x=x2;
}
    
template<class T,class A,class B,class C>
void Test_V_3(A& M,const B& b,const C& x)
{
    VECTOR<T,C::m> x0(x);
    ARRAY<T> x1(x),t2(x),t3(x);
    IDENTITY_ARRAY<> id(t2.m);
    auto x2=t2.Subset(id);
    ARRAY_VIEW<T> x3(t3);
    Test_V_4<T>(M,b,x0);
    Test_V_4<T>(M,b,x1);
    Test_V_4<T>(M,b,x2);
    Test_V_4<T>(M,b,x3);
}
    
template<class T,class A,class B,class C>
void Test_V_2(A& M,const B& b,const C& x)
{
    ARRAY<T> b1(b);
    auto b2=(b+b1)*(T).5;
    ARRAY_VIEW<T> b3(b1);
    Test_V_3<T>(M,b,x);
    Test_V_3<T>(M,b1,x);
    Test_V_3<T>(M,b2,x);
    Test_V_3<T>(M,b3,x);
}
    
template<class T,class A,class B,class C>
void Test_V_1(const A& M,const B& b,const C& x)
{
    A M0(M);
    MATRIX_MXN<T> M1(M);
    MATRIX_VIEW<T> M3(M1);
    Test_V_2<T>(M0,b,x);
    Test_V_2<T>(M1,b,x);
    Test_V_2<T>(M3,b,x);
}

template<class T,int m,int n>
void Test_V()
{
    RANDOM_NUMBERS<T> rand;

    MATRIX<T,m,n> M;
    VECTOR<T,m> b;
    VECTOR<T,n> x;
    rand.Fill_Uniform(M,-1,1);
    rand.Fill_Uniform(b,-1,1);
    rand.Fill_Uniform(x,-1,1);
    LOG::printf("orig: %P %P %P\n",M,b,x);
    Test_V_1<T>(M,b,x);
}

template<class T,class A,class B,class C>
void Test_M_4(A& M,const B& b,C& x)
{
    LOG::printf("try: %s %s %s\n",typeid(A).name(),typeid(B).name(),typeid(C).name());
    MATRIX_MXN<T> M2(M),x2(x);
    Least_Squares_Solve(M,b,x);
    LOG::printf("sol: %P   %P\n",(M2.Transpose_Times(M2*x-b)).Frobenius_Norm(),x);
    M=M2;
    x=x2;
    Least_Squares_Solve_Raw(M,b,x);
    LOG::printf("row: %P   %P\n",(M2.Transpose_Times(M2*x-b)).Frobenius_Norm(),x);
    M=M2;
    x=x2;
}
    
template<class T,class A,class B,class C>
void Test_M_3(A& M,const B& b,const C& x)
{
    C x0(x);
    MATRIX_MXN<T> x1(x);
    MATRIX_VIEW<T> x3(x1);
    Test_M_4<T>(M,b,x0);
    Test_M_4<T>(M,b,x1);
    Test_M_4<T>(M,b,x3);
}
    
template<class T,class A,class B,class C>
void Test_M_2(A& M,const B& b,const C& x)
{
    B b0(b);
    MATRIX_MXN<T> b1(b);
    MATRIX_VIEW<T> b3(b1);
    Test_M_3<T>(M,b0,x);
    Test_M_3<T>(M,b1,x);
    Test_M_3<T>(M,b3,x);
}
    
template<class T,class A,class B,class C>
void Test_M_1(const A& M,const B& b,const C& x)
{
    A M0(M);
    MATRIX_MXN<T> M1(M);
    MATRIX_VIEW<T> M3(M1);
    Test_M_2<T>(M0,b,x);
    Test_M_2<T>(M1,b,x);
    Test_M_2<T>(M3,b,x);
}

template<class T,int m,int n,int p>
void Test_M()
{
    RANDOM_NUMBERS<T> rand;

    MATRIX<T,m,n> M;
    MATRIX<T,m,p> b;
    MATRIX<T,n,p> x;
    rand.Fill_Uniform(M,-1,1);
    rand.Fill_Uniform(b,-1,1);
    rand.Fill_Uniform(x,-1,1);
    LOG::printf("orig: %P %P %P\n",M,b,x);
    Test_M_1<T>(M,b,x);
}

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    Test_V<double,1,1>();
    Test_V<double,1,2>();
    Test_V<double,1,3>();
    Test_V<double,1,4>();
    Test_V<double,2,1>();
    Test_V<double,2,2>();
    Test_V<double,2,3>();
    Test_V<double,2,4>();
    Test_V<double,3,1>();
    Test_V<double,3,2>();
    Test_V<double,3,3>();
    Test_V<double,3,4>();
    Test_V<double,4,1>();
    Test_V<double,4,2>();
    Test_V<double,4,3>();
    Test_V<double,4,4>();

    Test_M<double,2,2,2>();
    Test_M<double,2,3,2>();
    Test_M<double,3,2,2>();
    Test_M<double,3,3,2>();
    Test_M<double,2,2,3>();
    Test_M<double,2,3,3>();
    Test_M<double,3,2,3>();
    Test_M<double,3,3,3>();

    Test_V<float,1,1>();
    Test_V<float,1,2>();
    Test_V<float,1,3>();
    Test_V<float,1,4>();
    Test_V<float,2,1>();
    Test_V<float,2,2>();
    Test_V<float,2,3>();
    Test_V<float,2,4>();
    Test_V<float,3,1>();
    Test_V<float,3,2>();
    Test_V<float,3,3>();
    Test_V<float,3,4>();
    Test_V<float,4,1>();
    Test_V<float,4,2>();
    Test_V<float,4,3>();
    Test_V<float,4,4>();

    Test_M<float,2,2,2>();
    Test_M<float,2,3,2>();
    Test_M<float,3,2,2>();
    Test_M<float,3,3,2>();
    Test_M<float,2,2,3>();
    Test_M<float,2,3,3>();
    Test_M<float,3,2,3>();
    Test_M<float,3,3,3>();
    
    return 0;
}

