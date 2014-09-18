//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRIMITIVE_TENSORS
//##################################################################### 
#ifndef __PRIMITIVE_TENSORS__
#define __PRIMITIVE_TENSORS__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/PERMUTATION_TENSOR.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Tools/Tensors/TENSOR.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Tensors/VEC_ID_TENSOR_0.h>
#include <Tools/Tensors/VEC_ID_TENSOR_1.h>
#include <Tools/Tensors/VEC_ID_TENSOR_12.h>
#include <Tools/Tensors/VEC_ID_TENSOR_2.h>
#include <Tools/Tensors/ZERO_TENSOR.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

template<class T,int m,class TN>
ZERO_MATRIX<T,m> Contract_0(const TN& t,ZERO_VECTOR<T,m> z){return ZERO_MATRIX<T,m>();}

template<class T,int m>
MATRIX<T,m> Contract_0(const TENSOR<T,m>& t,const VECTOR<T,m>& v)
{
    MATRIX<T,m> M;
    for(int i=0;i<m;i++) M+=t.x(i)*v(i);
    return M;
}

template<class T,int m>
SYMMETRIC_MATRIX<T,m> Contract_0(const SYMMETRIC_TENSOR<T,m>& t,const VECTOR<T,m>& v)
{
    SYMMETRIC_MATRIX<T,m> M;
    for(int i=0;i<m;i++) M+=t.x(i)*v(i);
    return M;
}

template<class T,int m>
ZERO_MATRIX<T,m> Contract_0(const ZERO_TENSOR<T,m>& t,const VECTOR<T,m>& v) {return ZERO_MATRIX<T,m>();}

template<class T,int m>
SCALE_MATRIX<T,m> Contract_0(const VEC_ID_TENSOR_0<T,m>& t,const VECTOR<T,m>& v) {return SCALE_MATRIX<T,m>(t.v.Dot(v));}

template<class T,int m>
MATRIX<T,m> Contract_0(const VEC_ID_TENSOR_1<T,m>& t,const VECTOR<T,m>& v) {return MATRIX<T,m>::Outer_Product(t.v,v);}

template<class T,int m>
MATRIX<T,m> Contract_0(const VEC_ID_TENSOR_2<T,m>& t,const VECTOR<T,m>& v) {return MATRIX<T,m>::Outer_Product(v,t.v);}

template<class T,int m>
SYMMETRIC_MATRIX<T,m> Contract_0(const VEC_ID_TENSOR_12<T,m>& t,const VECTOR<T,m>& v) {return SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(t.v,v);}

template<class T,int m>
MATRIX<T,m> Contract_0(const PERMUTATION_TENSOR<T>& t,const VECTOR<T,m>& v) {return MATRIX<T,m>::Cross_Product_Matrix(-t.x*v);}

template<class T,int m,class MAT>
ZERO_TENSOR<T,m> Tensor_Product_0(const MAT& t,ZERO_VECTOR<T,m> z){return ZERO_TENSOR<T,m>();}

template<class T,int m>
TENSOR<T,m> Tensor_Product_0(const MATRIX<T,m>& M,const VECTOR<T,m>& v)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=M*v(i);
    return t;
}

template<class T,int m>
SYMMETRIC_TENSOR<T,m> Tensor_Product_0(const SYMMETRIC_MATRIX<T,m>& M,const VECTOR<T,m>& v)
{
    SYMMETRIC_TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=M*v(i);
    return t;
}

template<class T,int m>
ZERO_TENSOR<T,m> Tensor_Product_0(const ZERO_MATRIX<T,m>& M,const VECTOR<T,m>& v) {return ZERO_TENSOR<T,m>();}

template<class T,int m>
VEC_ID_TENSOR_0<T,m> Tensor_Product_0(IDENTITY_MATRIX<T,m> M,const VECTOR<T,m>& v) {return VEC_ID_TENSOR_0<T,m>(v);}

template<class T,int m>
VEC_ID_TENSOR_0<T,m> Tensor_Product_0(SCALE_MATRIX<T,m> M,const VECTOR<T,m>& v) {return VEC_ID_TENSOR_0<T,m>(M.x*v);}


template<class T,int m,class MAT>
ZERO_TENSOR<T,m> Tensor_Product_1(const MAT& t,ZERO_VECTOR<T,m> z){return ZERO_TENSOR<T,m>();}

template<class T,int m>
TENSOR<T,m> Tensor_Product_1(const MATRIX<T,m>& M,const VECTOR<T,m>& v)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=MATRIX<T,m>::Outer_Product(v,M.Row(i));
    return t;
}

template<class T,int m>
ZERO_TENSOR<T,m> Tensor_Product_1(const ZERO_MATRIX<T,m>& M,const VECTOR<T,m>& v) {return ZERO_TENSOR<T,m>();}

template<class T,int m>
VEC_ID_TENSOR_1<T,m> Tensor_Product_1(IDENTITY_MATRIX<T,m> M,const VECTOR<T,m>& v) {return VEC_ID_TENSOR_1<T,m>(v);}

template<class T,int m>
VEC_ID_TENSOR_1<T,m> Tensor_Product_1(SCALE_MATRIX<T,m> M,const VECTOR<T,m>& v) {return VEC_ID_TENSOR_1<T,m>(M.x*v);}


template<class T,int m,class MAT>
ZERO_TENSOR<T,m> Tensor_Product_2(const MAT& t,ZERO_VECTOR<T,m> z){return ZERO_TENSOR<T,m>();}

template<class T,int m>
TENSOR<T,m> Tensor_Product_2(const MATRIX<T,m>& M,const VECTOR<T,m>& v)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=MATRIX<T,m>::Outer_Product(M.Row(i),v);
    return t;
}

template<class T,int m>
ZERO_TENSOR<T,m> Tensor_Product_2(const ZERO_MATRIX<T,m>& M,const VECTOR<T,m>& v) {return ZERO_TENSOR<T,m>();}

template<class T,int m>
VEC_ID_TENSOR_2<T,m> Tensor_Product_2(IDENTITY_MATRIX<T,m> M,const VECTOR<T,m>& v) {return VEC_ID_TENSOR_2<T,m>(v);}

template<class T,int m>
VEC_ID_TENSOR_2<T,m> Tensor_Product_2(SCALE_MATRIX<T,m> M,const VECTOR<T,m>& v) {return VEC_ID_TENSOR_2<T,m>(M.x*v);}


template<class T,int m,class MAT>
ZERO_TENSOR<T,m> Symmetric_Tensor_Product_12(const MAT& t,ZERO_VECTOR<T,m> z){return ZERO_TENSOR<T,m>();}

template<class T,int m>
SYMMETRIC_TENSOR<T,m> Symmetric_Tensor_Product_12(const MATRIX<T,m>& M,const VECTOR<T,m>& v)
{
    SYMMETRIC_TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(M.Row(i),v);
    return t;
}

template<class T,int m>
ZERO_TENSOR<T,m> Symmetric_Tensor_Product_12(const ZERO_MATRIX<T,m>& M,const VECTOR<T,m>& v) {return ZERO_TENSOR<T,m>();}

template<class T,int m>
VEC_ID_TENSOR_12<T,m> Symmetric_Tensor_Product_12(IDENTITY_MATRIX<T,m> M,const VECTOR<T,m>& v) {return VEC_ID_TENSOR_12<T,m>(v);}

template<class T,int m>
VEC_ID_TENSOR_12<T,m> Symmetric_Tensor_Product_12(SCALE_MATRIX<T,m> M,const VECTOR<T,m>& v) {return VEC_ID_TENSOR_12<T,m>(M.x*v);}

template<class T,int m> TENSOR<T,m> operator+(const TENSOR<T,m>& a,const TENSOR<T,m>& b)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=a.x(i)+b.x(i);
    return t;
}
template<class T,int m> TENSOR<T,m> operator+(const TENSOR<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> TENSOR<T,m> operator+(const TENSOR<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b)
{
    TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) t.x(i)+=b.v(i);
    return t;
}
template<class T,int m> TENSOR<T,m> operator+(const TENSOR<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b)
{
    TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) for(int j=0;j<m;j++) t.x(i)(j,i)+=b.v(j);
    return t;
}
template<class T,int m> TENSOR<T,m> operator+(const TENSOR<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b)
{
    TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) for(int j=0;j<m;j++) t.x(i)(i,j)+=b.v(j);
    return t;
}
template<class T,int m> TENSOR<T,m> operator+(const TENSOR<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b)
{
    TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) for(int j=0;j<m;j++){t.x(i)(j,i)+=b.v(j);t.x(i)(i,j)+=b.v(j);}
    return t;
}
template<class T,int m> TENSOR<T,m> operator+(const TENSOR<T,m>& a,const PERMUTATION_TENSOR<T>& b)
{
    TENSOR<T,m> t=a;
    t.x(0)(1,2)+=b.x;
    t.x(1)(2,0)+=b.x;
    t.x(2)(0,1)+=b.x;
    t.x(0)(2,1)-=b.x;
    t.x(1)(0,2)-=b.x;
    t.x(2)(1,0)-=b.x;
    return t;
}

template<class T,int m> SYMMETRIC_TENSOR<T,m> operator+(const SYMMETRIC_TENSOR<T,m>& a,const SYMMETRIC_TENSOR<T,m>& b)
{
    SYMMETRIC_TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=a.x(i)+b.x(i);
    return t;
}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator+(const SYMMETRIC_TENSOR<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator+(const SYMMETRIC_TENSOR<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b)
{
    SYMMETRIC_TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) t.x(i)+=b.v(i);
    return t;
}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator+(const SYMMETRIC_TENSOR<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b)
{
    SYMMETRIC_TENSOR<T,m> t=a;
    for(int i=0;i<m;i++){
        for(int j=0;j<=i;j++) t.x(i).Element_Upper(j,i)+=b.v(j);
        for(int j=i;j<m;j++) t.x(i).Element_Lower(j,i)+=b.v(j);}
    return t;
}

template<class T,int m> TENSOR<T,m> operator+(const ZERO_TENSOR<T,m>& a,const TENSOR<T,m>& b){return b;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator+(const ZERO_TENSOR<T,m>& a,const SYMMETRIC_TENSOR<T,m>& b){return b;}
template<class T,int m> ZERO_TENSOR<T,m> operator+(const ZERO_TENSOR<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> VEC_ID_TENSOR_0<T,m> operator+(const ZERO_TENSOR<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return b;}
template<class T,int m> VEC_ID_TENSOR_1<T,m> operator+(const ZERO_TENSOR<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return b;}
template<class T,int m> VEC_ID_TENSOR_2<T,m> operator+(const ZERO_TENSOR<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return b;}
template<class T,int m> VEC_ID_TENSOR_12<T,m> operator+(const ZERO_TENSOR<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return b;}
template<class T,int m> PERMUTATION_TENSOR<T> operator+(const ZERO_TENSOR<T,m>& a,const PERMUTATION_TENSOR<T>& b){return b;}

template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_0<T,m>& a,const TENSOR<T,m>& b){return b+a;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator+(const VEC_ID_TENSOR_0<T,m>& a,const SYMMETRIC_TENSOR<T,m>& b){return b+a;}
template<class T,int m> VEC_ID_TENSOR_0<T,m> operator+(const VEC_ID_TENSOR_0<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> VEC_ID_TENSOR_0<T,m> operator+(const VEC_ID_TENSOR_0<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return VEC_ID_TENSOR_0<T,m>(a.v+b.v);}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_0<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return TENSOR<T,m>()+a+b;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_0<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return TENSOR<T,m>()+a+b;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator+(const VEC_ID_TENSOR_0<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return SYMMETRIC_TENSOR<T,m>()+a+b;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_0<T,m>& a,const PERMUTATION_TENSOR<T>& b){return TENSOR<T,m>()+a+b;}

template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_1<T,m>& a,const TENSOR<T,m>& b){return b+a;}
template<class T,int m> VEC_ID_TENSOR_1<T,m> operator+(const VEC_ID_TENSOR_1<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_1<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return b+a;}
template<class T,int m> VEC_ID_TENSOR_1<T,m> operator+(const VEC_ID_TENSOR_1<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return VEC_ID_TENSOR_1<T,m>(a.v+b.v);}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_1<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return TENSOR<T,m>()+a+b;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_1<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return TENSOR<T,m>()+a+b;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_1<T,m>& a,const PERMUTATION_TENSOR<T>& b){return TENSOR<T,m>()+a+b;}

template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_2<T,m>& a,const TENSOR<T,m>& b){return b+a;}
template<class T,int m> VEC_ID_TENSOR_2<T,m> operator+(const VEC_ID_TENSOR_2<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_2<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return b+a;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_2<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return b+a;}
template<class T,int m> VEC_ID_TENSOR_2<T,m> operator+(const VEC_ID_TENSOR_2<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return VEC_ID_TENSOR_2<T,m>(a.v+b.v);}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_2<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return TENSOR<T,m>()+a+b;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_2<T,m>& a,const PERMUTATION_TENSOR<T>& b){return TENSOR<T,m>()+a+b;}

template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_12<T,m>& a,const TENSOR<T,m>& b){return b+a;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator+(const VEC_ID_TENSOR_12<T,m>& a,const SYMMETRIC_TENSOR<T,m>& b){return b+a;}
template<class T,int m> VEC_ID_TENSOR_12<T,m> operator+(const VEC_ID_TENSOR_12<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator+(const VEC_ID_TENSOR_12<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return b+a;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_12<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return b+a;}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_12<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return b+a;}
template<class T,int m> VEC_ID_TENSOR_12<T,m> operator+(const VEC_ID_TENSOR_12<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return VEC_ID_TENSOR_12<T,m>(a.v+b.v);}
template<class T,int m> TENSOR<T,m> operator+(const VEC_ID_TENSOR_12<T,m>& a,const PERMUTATION_TENSOR<T>& b){return TENSOR<T,m>()+a+b;}

template<class T,int m> TENSOR<T,m> operator+(const PERMUTATION_TENSOR<T>& a,const TENSOR<T,m>& b){return b+a;}
template<class T,int m> PERMUTATION_TENSOR<T> operator+(const PERMUTATION_TENSOR<T>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> TENSOR<T,m> operator+(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR_0<T,m>& b){return b+a;}
template<class T,int m> TENSOR<T,m> operator+(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR_1<T,m>& b){return b+a;}
template<class T,int m> TENSOR<T,m> operator+(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR_2<T,m>& b){return b+a;}
template<class T,int m> TENSOR<T,m> operator+(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR_12<T,m>& b){return b+a;}
template<class T,int m> PERMUTATION_TENSOR<T> operator+(const PERMUTATION_TENSOR<T>& a,const PERMUTATION_TENSOR<T>& b){return PERMUTATION_TENSOR<T>(a.x+b.x);}

template<class T,int m> TENSOR<T,m> operator-(const TENSOR<T,m>& a,const TENSOR<T,m>& b)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=a.x(i)-b.x(i);
    return t;
}
template<class T,int m> TENSOR<T,m> operator-(const TENSOR<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> TENSOR<T,m> operator-(const TENSOR<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b)
{
    TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) t.x(i)-=b.v(i);
    return t;
}
template<class T,int m> TENSOR<T,m> operator-(const TENSOR<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b)
{
    TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) for(int j=0;j<m;j++) t.x(i)(j,i)-=b.v(j);
    return t;
}
template<class T,int m> TENSOR<T,m> operator-(const TENSOR<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b)
{
    TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) for(int j=0;j<m;j++) t.x(i)(i,j)-=b.v(j);
    return t;
}
template<class T,int m> TENSOR<T,m> operator-(const TENSOR<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b)
{
    TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) for(int j=0;j<m;j++){t.x(i)(j,i)-=b.v(j);t.x(i)(i,j)-=b.v(j);}
    return t;
}
template<class T,int m> TENSOR<T,m> operator-(const TENSOR<T,m>& a,const PERMUTATION_TENSOR<T>& b)
{
    TENSOR<T,m> t=a;
    t.x(0)(1,2)-=b.x;
    t.x(1)(2,0)-=b.x;
    t.x(2)(0,1)-=b.x;
    t.x(0)(2,1)+=b.x;
    t.x(1)(0,2)+=b.x;
    t.x(2)(1,0)+=b.x;
    return t;
}

template<class T,int m> SYMMETRIC_TENSOR<T,m> operator-(const SYMMETRIC_TENSOR<T,m>& a,const SYMMETRIC_TENSOR<T,m>& b)
{
    SYMMETRIC_TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=a.x(i)-b.x(i);
    return t;
}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator-(const SYMMETRIC_TENSOR<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator-(const SYMMETRIC_TENSOR<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b)
{
    SYMMETRIC_TENSOR<T,m> t=a;
    for(int i=0;i<m;i++) t.x(i)-=b.v(i);
    return t;
}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator-(const SYMMETRIC_TENSOR<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b)
{
    SYMMETRIC_TENSOR<T,m> t=a;
    for(int i=0;i<m;i++){
        for(int j=0;j<=i;j++) t.x(i).Element_Upper(j,i)-=b.v(j);
        for(int j=i;j<m;j++) t.x(i).Element_Lower(j,i)-=b.v(j);}
    return t;
}

template<class T,int m> TENSOR<T,m> operator-(const ZERO_TENSOR<T,m>& a,const TENSOR<T,m>& b){return -b;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator-(const ZERO_TENSOR<T,m>& a,const SYMMETRIC_TENSOR<T,m>& b){return -b;}
template<class T,int m> ZERO_TENSOR<T,m> operator-(const ZERO_TENSOR<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> VEC_ID_TENSOR_0<T,m> operator-(const ZERO_TENSOR<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return -b;}
template<class T,int m> VEC_ID_TENSOR_1<T,m> operator-(const ZERO_TENSOR<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return -b;}
template<class T,int m> VEC_ID_TENSOR_2<T,m> operator-(const ZERO_TENSOR<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return -b;}
template<class T,int m> VEC_ID_TENSOR_12<T,m> operator-(const ZERO_TENSOR<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return -b;}
template<class T,int m> PERMUTATION_TENSOR<T> operator-(const ZERO_TENSOR<T,m>& a,const PERMUTATION_TENSOR<T>& b){return -b;}

template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_0<T,m>& a,const TENSOR<T,m>& b){return -b+a;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator-(const VEC_ID_TENSOR_0<T,m>& a,const SYMMETRIC_TENSOR<T,m>& b){return -b+a;}
template<class T,int m> VEC_ID_TENSOR_0<T,m> operator-(const VEC_ID_TENSOR_0<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> VEC_ID_TENSOR_0<T,m> operator-(const VEC_ID_TENSOR_0<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return VEC_ID_TENSOR_0<T,m>(a.v-b.v);}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_0<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_0<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator-(const VEC_ID_TENSOR_0<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return SYMMETRIC_TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_0<T,m>& a,const PERMUTATION_TENSOR<T>& b){return TENSOR<T,m>()+a-b;}

template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_1<T,m>& a,const TENSOR<T,m>& b){return -b+a;}
template<class T,int m> VEC_ID_TENSOR_1<T,m> operator-(const VEC_ID_TENSOR_1<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_1<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> VEC_ID_TENSOR_1<T,m> operator-(const VEC_ID_TENSOR_1<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return VEC_ID_TENSOR_1<T,m>(a.v-b.v);}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_1<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_1<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_1<T,m>& a,const PERMUTATION_TENSOR<T>& b){return TENSOR<T,m>()+a-b;}

template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_2<T,m>& a,const TENSOR<T,m>& b){return -b+a;}
template<class T,int m> VEC_ID_TENSOR_2<T,m> operator-(const VEC_ID_TENSOR_2<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_2<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_2<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> VEC_ID_TENSOR_2<T,m> operator-(const VEC_ID_TENSOR_2<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return VEC_ID_TENSOR_2<T,m>(a.v-b.v);}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_2<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_2<T,m>& a,const PERMUTATION_TENSOR<T>& b){return TENSOR<T,m>()+a-b;}

template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_12<T,m>& a,const TENSOR<T,m>& b){return -b+a;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator-(const VEC_ID_TENSOR_12<T,m>& a,const SYMMETRIC_TENSOR<T,m>& b){return -b+a;}
template<class T,int m> VEC_ID_TENSOR_12<T,m> operator-(const VEC_ID_TENSOR_12<T,m>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> SYMMETRIC_TENSOR<T,m> operator-(const VEC_ID_TENSOR_12<T,m>& a,const VEC_ID_TENSOR_0<T,m>& b){return SYMMETRIC_TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_12<T,m>& a,const VEC_ID_TENSOR_1<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_12<T,m>& a,const VEC_ID_TENSOR_2<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> VEC_ID_TENSOR_12<T,m> operator-(const VEC_ID_TENSOR_12<T,m>& a,const VEC_ID_TENSOR_12<T,m>& b){return VEC_ID_TENSOR_12<T,m>(a.v-b.v);}
template<class T,int m> TENSOR<T,m> operator-(const VEC_ID_TENSOR_12<T,m>& a,const PERMUTATION_TENSOR<T>& b){return TENSOR<T,m>()+a-b;}

template<class T,int m> TENSOR<T,m> operator-(const PERMUTATION_TENSOR<T>& a,const TENSOR<T,m>& b){return -b+a;}
template<class T,int m> PERMUTATION_TENSOR<T> operator-(const PERMUTATION_TENSOR<T>& a,const ZERO_TENSOR<T,m>& b){return a;}
template<class T,int m> TENSOR<T,m> operator-(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR_0<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR_1<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR_2<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> TENSOR<T,m> operator-(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR_12<T,m>& b){return TENSOR<T,m>()+a-b;}
template<class T,int m> PERMUTATION_TENSOR<T> operator-(const PERMUTATION_TENSOR<T>& a,const PERMUTATION_TENSOR<T>& b){return PERMUTATION_TENSOR<T>(a.x-b.x);}

template<class T,int m> SYMMETRIC_TENSOR<T,m> Contract_0(const SYMMETRIC_TENSOR<T,m>& a,const MATRIX<T,m>& M)
{
    SYMMETRIC_TENSOR<T,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            t.x(j)+=a.x(i)*M(i,j);
    return t;
}
template<class T,int m> TENSOR<T,m> Contract_0(const TENSOR<T,m>& a,const MATRIX<T,m>& M)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            t.x(j)+=a.x(i)*M(i,j);
    return t;
}
template<class T,int m> ZERO_TENSOR<T,m> Contract_0(const ZERO_TENSOR<T,m>& a,const MATRIX<T,m>& M){return a;}
template<class T,int m> VEC_ID_TENSOR_0<T,m> Contract_0(const VEC_ID_TENSOR_0<T,m>& a,const MATRIX<T,m>& M){return VEC_ID_TENSOR_0<T,m>(M.Transpose_Times(a.v));}
template<class T,int m> TENSOR<T,m> Contract_0(const VEC_ID_TENSOR_1<T,m>& a,const MATRIX<T,m>& M){return Tensor_Product_1(M.Transposed(),a.v);}
template<class T,int m> TENSOR<T,m> Contract_0(const VEC_ID_TENSOR_2<T,m>& a,const MATRIX<T,m>& M){return Tensor_Product_2(M.Transposed(),a.v);}
template<class T,int m> SYMMETRIC_TENSOR<T,m> Contract_0(const VEC_ID_TENSOR_12<T,m>& a,const MATRIX<T,m>& M){return Symmetric_Tensor_Product_12(M.Transposed(),a.v);}
template<class T,int m> TENSOR<T,m> Contract_0(const PERMUTATION_TENSOR<T>& a,const MATRIX<T,m>& M)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=MATRIX<T,m>::Cross_Product_Matrix(-a.x*M.Column(i));
    return t;
}

template<class T,int m> ZERO_TENSOR<T,m> Contract_1(const ZERO_TENSOR<T,m>& p,const MATRIX<T,m>& M) {return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Contract_1(const ZERO_TENSOR<T,m>& p,const ZERO_MATRIX<T,m>& M) {return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Contract_1(const ZERO_TENSOR<T,m>& p,IDENTITY_MATRIX<T,m> M) {return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Contract_1(const ZERO_TENSOR<T,m>& p,SCALE_MATRIX<T,m> M) {return ZERO_TENSOR<T,m>();}

template<class T,int m> TENSOR<T,m> Contract_1(const TENSOR<T,m>& p,const MATRIX<T,m>& M)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=M.Transpose_Times(p.x(i));
    return t;
}
template<class T,int m> ZERO_TENSOR<T,m> Contract_1(const TENSOR<T,m>& p,const ZERO_MATRIX<T,m>& M) {return ZERO_TENSOR<T,m>();}
template<class T,int m> TENSOR<T,m> Contract_1(const TENSOR<T,m>& p,IDENTITY_MATRIX<T,m> M) {return p;}
template<class T,int m> TENSOR<T,m> Contract_1(const TENSOR<T,m>& p,SCALE_MATRIX<T,m> M) {return p*M.x;}

template<class T,int m> TENSOR<T,m> Contract_1(const PERMUTATION_TENSOR<T>& p,const MATRIX<T,m>& M)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++){
        MATRIX<T,m> c=MATRIX<T,m>::Cross_Product_Matrix(p.x*M.Column(i));
        for(int j=0;j<m;j++) t.x(j).Set_Row(i,c.Row(j));}
    return t;
}
template<class T,int m> ZERO_TENSOR<T,m> Contract_1(const PERMUTATION_TENSOR<T>& p,const ZERO_MATRIX<T,m>& M) {return ZERO_TENSOR<T,m>();}
template<class T,int m> PERMUTATION_TENSOR<T> Contract_1(const PERMUTATION_TENSOR<T>& p,IDENTITY_MATRIX<T,m> M) {return p;}
template<class T,int m> PERMUTATION_TENSOR<T> Contract_1(const PERMUTATION_TENSOR<T>& p,SCALE_MATRIX<T,m> M) {return p*M.x;}

template<class T,int m> TENSOR<T,m> Contract_2(const TENSOR<T,m>& p,const MATRIX<T,m>& M)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++) t.x(i)=p.x(i)*m;
    return t;
}
template<class T,int m> TENSOR<T,m> Contract_2(const PERMUTATION_TENSOR<T>& p,const MATRIX<T,m>& M)
{
    TENSOR<T,m> t;
    for(int i=0;i<m;i++){
        MATRIX<T,m> c=MATRIX<T,m>::Cross_Product_Matrix(p.x*M.Column(i));
        for(int j=0;j<m;j++) t.x(j).Set_Column(i,c.Column(j));}
    return t;
}
template<class T,int m> ZERO_TENSOR<T,m> Contract_2(const PERMUTATION_TENSOR<T>& p,const ZERO_MATRIX<T,m>& M) {return ZERO_TENSOR<T,m>();}
template<class T,int m> PERMUTATION_TENSOR<T> Contract_2(const PERMUTATION_TENSOR<T>& p,IDENTITY_MATRIX<T,m> M) {return p;}
template<class T,int m> PERMUTATION_TENSOR<T> Contract_2(const PERMUTATION_TENSOR<T>& p,SCALE_MATRIX<T,m> M) {return p*M.x;}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const ZERO_MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const ZERO_MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const ZERO_MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const ZERO_MATRIX<T,m>& cm1,const MATRIX<T,m>){return ZERO_TENSOR<T,m>();}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const IDENTITY_MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const IDENTITY_MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const IDENTITY_MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> SYMMETRIC_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const IDENTITY_MATRIX<T,m>& cm1,const MATRIX<T,m>& cm2)
{return Contract_2(t,cm2).Twice_Symmetric_Part_12();}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const SCALE_MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const SCALE_MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const SCALE_MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> SYMMETRIC_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const SCALE_MATRIX<T,m>& cm1,const MATRIX<T,m>& cm2)
{return Symmetric_Double_Contract_12(t*cm1.x,IDENTITY_MATRIX<T,m>(),cm2);}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> SYMMETRIC_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2)
{return Symmetric_Double_Contract_12(-t,cm2,cm1);}
template<class T,int m> SYMMETRIC_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2)
{return Symmetric_Double_Contract_12(-t,cm2,cm1);}
template<class T,int m> SYMMETRIC_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const MATRIX<T,m>& cm1,const MATRIX<T,m>& cm2)
{
    SYMMETRIC_TENSOR<T,m> s;
    VECTOR<T,m> c10=cm1.Row(0),c11=cm1.Row(1),c12=cm1.Row(2),c20=t.x*cm2.Row(0),c21=t.x*cm2.Row(1),c22=t.x*cm2.Row(2);
    s.x(0)=SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c11,c22)-SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c12,c21);
    s.x(1)=SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c12,c20)-SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c10,c22);
    s.x(2)=SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c10,c21)-SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c11,c20);
    return s;
}

template<class T_TEN> typename ENABLE_IF<IS_TENSOR<T_TEN>::value,T_TEN>::TYPE Choose(const T_TEN& a,const T_TEN& b);
template<class T_TEN,class T_TEN1> typename ENABLE_IF<IS_TENSOR<T_TEN>::value&&IS_TENSOR<T_TEN1>::value&&(!IS_SYM_TENSOR<T_TEN>::value || !IS_SYM_TENSOR<T_TEN1>::value),TENSOR<typename T_TEN::SCALAR,T_TEN::m> >::TYPE Choose(const T_TEN& a,const T_TEN1& b);
template<class T_TEN,class T_TEN1> typename ENABLE_IF<IS_SYM_TENSOR<T_TEN>::value && IS_SYM_TENSOR<T_TEN1>::value,SYMMETRIC_TENSOR<typename T_TEN::SCALAR,T_TEN::m> >::TYPE Choose(const T_TEN& a,const T_TEN1& b);
template<class T,int m,class T_TEN> T_TEN Choose(const T_TEN& a,const ZERO_TENSOR<T,m>& b);
template<class T,int m,class T_TEN> T_TEN Choose(const ZERO_TENSOR<T,m>& a,const T_TEN& b);
template<class T,int m> ZERO_TENSOR<T,m> Choose(const ZERO_TENSOR<T,m>& a,const ZERO_TENSOR<T,m>& b);

template<class T_TEN> typename ENABLE_IF<IS_TENSOR<T_TEN>::value>::TYPE Fill_From(T_TEN& a,const T_TEN& b){a=b;}
template<class T_TEN,class T_TEN1> typename ENABLE_IF<IS_TENSOR<T_TEN>::value&&IS_TENSOR<T_TEN1>::value>::TYPE Fill_From(T_TEN& a,const T_TEN1& b){a=b;}
}
#endif
