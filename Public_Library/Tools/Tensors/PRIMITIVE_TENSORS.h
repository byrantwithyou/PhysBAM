//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRIMITIVE_TENSORS
//##################################################################### 
#ifndef __PRIMITIVE_TENSORS__
#define __PRIMITIVE_TENSORS__

#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Tensors/DIAGONAL_TENSOR.h>
#include <Tools/Tensors/PERMUTATION_TENSOR.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Tools/Tensors/TENSOR.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Tensors/VEC_ID_SYM_TENSOR.h>
#include <Tools/Tensors/VEC_ID_TENSOR.h>
#include <Tools/Tensors/ZERO_TENSOR.h>
#include <cmath>
namespace PhysBAM{

/*
ZERO_VECTOR
VECTOR

ZERO_MATRIX
IDENTITY_MATRIX
SCALE_MATRIX
SYMMETRIC_MATRIX
MATRIX

ZERO_TENSOR
TENSOR
SYMMETRIC_TENSOR
VEC_ID_TENSOR
VEC_ID_SYM_TENSOR
PERMUTATION_TENSOR
*/

//#####################################################################
// Function Contract<s>(T,V)
//##################################################################### 

template<int s,class T,class TN>
typename enable_if<IS_TENSOR<TN>::value,ZERO_MATRIX<T,s==0?TN::n:TN::m,s==2?TN::n:TN::p> >::type
Contract(const TN& t,ZERO_VECTOR<T,(s==0?TN::m:s==1?TN::n:TN::p)> z)
{return ZERO_MATRIX<T,s==0?TN::n:TN::m,s==2?TN::n:TN::p>();}

template<int s,class T,int m,int n,int p> ZERO_MATRIX<T,s==0?n:m,s==2?n:p>
Contract(const ZERO_TENSOR<T,m,n,p>& t,const VECTOR<T,(s==0?m:s==1?n:p)>& v)
{return ZERO_MATRIX<T,s==0?n:m,s==2?n:p>();}

template<int s,class T,int m,int n,int p>
typename enable_if<s==0,MATRIX<T,n,p> >::type
Contract(const TENSOR<T,m,n,p>& t,const VECTOR<T,m>& v)
{
    MATRIX<T,n,p> M;
    for(int i=0;i<m;i++) M+=t.x(i)*v(i);
    return M;
}

template<int s,class T,int m,int n,int p>
typename enable_if<s==1,MATRIX<T,m,p> >::type
Contract(const TENSOR<T,m,n,p>& t,const VECTOR<T,n>& v)
{
    MATRIX<T,m,p> M;
    for(int i=0;i<m;i++) M.Set_Row(i,t.x(i).Transpose_Times(v));
    return M;
}

template<int s,class T,int m,int n,int p>
typename enable_if<s==2,MATRIX<T,m,n> >::type
Contract(const TENSOR<T,m,n,p>& t,const VECTOR<T,p>& v)
{
    MATRIX<T,m,n> M;
    for(int i=0;i<m;i++) M.Set_Row(i,t.x(i)*v);
    return M;
}

template<int s,class T,int m,int n> SYMMETRIC_MATRIX<T,n>
Contract(const SYMMETRIC_TENSOR<T,s,m,n>& t,const VECTOR<T,m>& v)
{
    SYMMETRIC_MATRIX<T,n> M;
    for(int i=0;i<m;i++) M+=t.x(i)*v(i);
    return M;
}

template<int s,class T,int m,int n> MATRIX<T,m,n>
Contract(const SYMMETRIC_TENSOR<T,s==0?1:0,m,n>& t,const VECTOR<T,n>& v)
{
    MATRIX<T,m,n> M;
    for(int i=0;i<m;i++) M.Set_Row(i,t.x(i).Transpose_Times(v));
    return M;
}

template<int s,class T,int m,int n> MATRIX<T,n,m>
Contract(const SYMMETRIC_TENSOR<T,s==2?1:2,m,n>& t,const VECTOR<T,n>& v)
{
    MATRIX<T,n,m> M;
    for(int i=0;i<m;i++) M.Set_Column(i,t.x(i)*v);
    return M;
}

template<int s,class T,int m,int n> SCALE_MATRIX<T,n>
Contract(const VEC_ID_TENSOR<T,s,m,n>& t,const VECTOR<T,m>& v)
{
    return SCALE_MATRIX<T,n>(t.v.Dot(v));
}

template<int s,class T,int m,int n> MATRIX<T,m,n>
Contract(const VEC_ID_TENSOR<T,s==0?1:0,m,n>& t,const VECTOR<T,n>& v)
{
    return MATRIX<T,m,n>::Outer_Product(t.v,v);
}

template<int s,class T,int m,int n> MATRIX<T,n,m>
Contract(const VEC_ID_TENSOR<T,s==2?1:2,m,n>& t,const VECTOR<T,n>& v)
{
    return MATRIX<T,n,m>::Outer_Product(v,t.v);
}

template<int s,class T,int m> SYMMETRIC_MATRIX<T,m>
Contract(const VEC_ID_SYM_TENSOR<T,s,m>& t,const VECTOR<T,m>& v)
{
    return SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(v,t.v);
}

template<int s,class T,int m> MATRIX<T,m>
Contract(const VEC_ID_SYM_TENSOR<T,s==0?1:0,m>& t,const VECTOR<T,m>& v)
{
    return MATRIX<T,m>::Outer_Product(v,t.v)+v.Dot(t.v);
}

template<int s,class T,int m> MATRIX<T,m>
Contract(const VEC_ID_SYM_TENSOR<T,s==2?1:2,m>& t,const VECTOR<T,m>& v)
{
    return MATRIX<T,m>::Outer_Product(t.v,v)+v.Dot(t.v);
}

template<int s,class T> MATRIX<T,3>
Contract(const PERMUTATION_TENSOR<T>& t,const VECTOR<T,3>& v)
{
    return MATRIX<T,3>::Cross_Product_Matrix(s==1?t.x*v:-t.x*v);
}

template<int s,class T,int m> DIAGONAL_MATRIX<T,m>
Contract(const DIAGONAL_TENSOR<T,m>& t,const VECTOR<T,m>& v)
{
    return DIAGONAL_MATRIX<T,m>(v*t.v);
}

//#####################################################################
// Function Contract<r,s>(T)
//##################################################################### 

template<int r,int s,class T,int m,int n,int p> typename enable_if<(r+s!=1 || m==n)&&(r+s!=2 || m==p)&&(r+s!=3 || n==p),ZERO_VECTOR<T,r+s==1?p:r+s==2?n:m> >::type
Contract(const ZERO_TENSOR<T,m,n,p>& t)
{return ZERO_VECTOR<T,r+s==1?p:r+s==2?n:m>();}

template<int r,int s,class T,int m,int p>
typename enable_if<r+s==1,VECTOR<T,p> >::type
Contract(const TENSOR<T,m,m,p>& t)
{
    VECTOR<T,p> v;
    for(int i=0;i<m;i++) v+=t.x(i).Row(i);
    return v;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r+s==2,VECTOR<T,n> >::type
Contract(const TENSOR<T,m,n,m>& t)
{
    VECTOR<T,n> v;
    for(int i=0;i<m;i++) v+=t.x(i).Column(i);
    return v;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r+s==3,VECTOR<T,m> >::type
Contract(const TENSOR<T,m,n,n>& t)
{
    VECTOR<T,m> v;
    for(int i=0;i<m;i++) v(i)=t.x(i).Trace();
    return v;
}

template<int r,int s,class T,int m,int n> VECTOR<T,m>
Contract(const SYMMETRIC_TENSOR<T,3-r-s,m,n>& t)
{
    VECTOR<T,m> v;
    for(int i=0;i<m;i++) v(i)=t.x(i).Trace();
    return v;
}

template<int r,int s,class T,int m,int u>
typename enable_if<r+s+u!=3,VECTOR<T,m> >::type
Contract(const SYMMETRIC_TENSOR<T,u,m>& t)
{
    VECTOR<T,m> v;
    for(int i=0;i<m;i++) v+=t.x(i).Column(i);
    return v;
}

template<int r,int s,class T,int m,int n> VECTOR<T,m>
Contract(const VEC_ID_TENSOR<T,3-r-s,m,n>& t)
{
    return t.v*n;
}

template<int r,int s,class T,int m,int u>
typename enable_if<r+s+u!=3,VECTOR<T,m> >::type
Contract(const VEC_ID_TENSOR<T,u,m,m>& t)
{
    return t.v;
}

template<int r,int s,class T,int m> VECTOR<T,m>
Contract(const VEC_ID_SYM_TENSOR<T,3-r-s,m>& t)
{
    return t.v*2;
}

template<int r,int s,class T,int u,int m>
typename enable_if<u!=3-r-s,VECTOR<T,m> >::type
Contract(const VEC_ID_SYM_TENSOR<T,u,m>& t)
{
    return t.v*(m+1);
}

template<int r,int s,class T> ZERO_VECTOR<T,3>
Contract(const PERMUTATION_TENSOR<T>& t)
{
    return ZERO_VECTOR<T,3>();
}

//#####################################################################
// Function Tensor_Product<s>(M,V)
//##################################################################### 

template<int s,class MAT,class T,int m>
typename enable_if<IS_MATRIX<MAT>::value,ZERO_TENSOR<T,s==0?m:MAT::m,s==1?m:s==2?MAT::n:MAT::m,s==2?m:MAT::n> >::type
Tensor_Product(const MAT&,const ZERO_VECTOR<T,m>&)
{
    return ZERO_TENSOR<T,s==0?m:MAT::m,s==1?m:s==2?MAT::n:MAT::m,s==2?m:MAT::n>();
}

template<int s,class T,int m,int n,int p>
ZERO_TENSOR<T,s==0?p:m,s==1?p:s==2?n:m,s==2?p:n>
Tensor_Product(const ZERO_MATRIX<T,m,n>&,const VECTOR<T,p>&)
{
    return ZERO_TENSOR<T,s==0?p:m,s==1?p:s==2?n:m,s==2?p:n>();
}

template<int s,class T,int m,int n,int p>
typename enable_if<s==0,TENSOR<T,p,m,n> >::type
Tensor_Product(const MATRIX<T,m,n>& M,const VECTOR<T,p>& v)
{
    TENSOR<T,p,m,n> t;
    for(int i=0;i<p;i++) t.x(i)=M*v(i);
    return t;
}

template<int s,class T,int m,int n,int p>
typename enable_if<s==1,TENSOR<T,m,p,n> >::type
Tensor_Product(const MATRIX<T,m,n>& M,const VECTOR<T,p>& v)
{
    TENSOR<T,m,p,n> t;
    for(int i=0;i<m;i++) t.x(i)=MATRIX<T,p,n>::Outer_Product(v,M.Row(i));
    return t;
}

template<int s,class T,int m,int n,int p>
typename enable_if<s==2,TENSOR<T,m,n,p> >::type
Tensor_Product(const MATRIX<T,m,n>& M,const VECTOR<T,p>& v)
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<m;i++) t.x(i)=MATRIX<T,n,p>::Outer_Product(M.Row(i),v);
    return t;
}

template<int s,class T,int m,int n>
SYMMETRIC_TENSOR<T,s,n,m>
Tensor_Product(const SYMMETRIC_MATRIX<T,m>& M,const VECTOR<T,n>& v)
{
    SYMMETRIC_TENSOR<T,s,n,m> t;
    for(int i=0;i<n;i++) t.x(i)=M*v(i);
    return t;
}

template<int s,class T,int m,int n>
SYMMETRIC_TENSOR<T,s,n,m>
Tensor_Product(const DIAGONAL_MATRIX<T,m>& M,const VECTOR<T,n>& v)
{
    SYMMETRIC_TENSOR<T,s,n,m> t;
    for(int i=0;i<n;i++) t.x(i)=M*v(i);
    return t;
}

template<int s,class T,int m,int n>
VEC_ID_TENSOR<T,s,n,m>
Tensor_Product(const SCALE_MATRIX<T,m>& M,const VECTOR<T,n>& v)
{
    return VEC_ID_TENSOR<T,s,n,m>(M.x*v);
}

template<int s,class T,int m,int n>
VEC_ID_TENSOR<T,s,n,m>
Tensor_Product(const IDENTITY_MATRIX<T,m>& M,const VECTOR<T,n>& v)
{
    return VEC_ID_TENSOR<T,s,n,m>(v);
}

//#####################################################################
// Function Transposed<r,s>(T)
//##################################################################### 

template<int r,int s,class T,int m,int n,int p>
ZERO_TENSOR<T,r+s==1?n:r+s==2?p:m,r+s==1?m:r+s==2?n:p,r+s==1?p:r+s==2?m:n>
Transposed(const ZERO_TENSOR<T,m,n,p>& t)
{
    return ZERO_TENSOR<T,r+s==1?n:r+s==2?p:m,r+s==1?m:r+s==2?n:p,r+s==1?p:r+s==2?m:n>();
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<r+s==1,TENSOR<T,n,m,p> >::type
Transposed(const TENSOR<T,m,n,p>& t)
{
    TENSOR<T,n,m,p> z;
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) z.x(j).Set_Row(i,t.x(i).Row(j));
    return z;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<r+s==2,TENSOR<T,p,n,m> >::type
Transposed(const TENSOR<T,m,n,p>& t)
{
    TENSOR<T,p,n,m> z;
    for(int i=0;i<m;i++) for(int j=0;j<p;j++) z.x(j).Set_Column(i,t.x(i).Column(j));
    return z;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<r+s==3,TENSOR<T,m,p,n> >::type
Transposed(const TENSOR<T,m,n,p>& t)
{
    TENSOR<T,m,p,n> z;
    for(int i=0;i<m;i++) z.x(i)=t.x(i).Transposed();
    return z;
}

template<int r,int s,class TEN>
typename enable_if<IS_SYM_TENSOR<3-r-s,TEN>::value,TEN>::type
Transposed(const TEN& t)
{
    return t;
}

template<int r,int s,class T,int m,int n,int u>
typename enable_if<r+s+u!=3,SYMMETRIC_TENSOR<T,r+s-u,m,n> >::type
Transposed(const SYMMETRIC_TENSOR<T,u,m,n>& t)
{
    SYMMETRIC_TENSOR<T,r+s-u,m,n> z;
    z.x=t.x;
    return z;
}

template<int r,int s,class T,int m,int n,int u>
typename enable_if<r+s+u!=3,VEC_ID_TENSOR<T,r+s-u,m,n> >::type
Transposed(const VEC_ID_TENSOR<T,u,m,n>& t)
{
    return VEC_ID_TENSOR<T,r+s-u,m,n>(t.v);
}

template<int r,int s,class T,int u,int m>
typename enable_if<r+s+u!=3,VEC_ID_SYM_TENSOR<T,r+s-u,m> >::type
Transposed(const VEC_ID_SYM_TENSOR<T,u,m>& t)
{
    return VEC_ID_SYM_TENSOR<T,r+s-u,m>(t.v);
}

template<int r,int s,class T> PERMUTATION_TENSOR<T>
Transposed(const PERMUTATION_TENSOR<T>& t)
{
    return -t;
}

//#####################################################################
// Function Twice_Symmetric_Part<r,s>(T)
//##################################################################### 

template<int r,int s,class TEN>
typename enable_if<IS_SYM_TENSOR<3-r-s,TEN>::value,decltype(TEN()*2)>::type
Twice_Symmetric_Part(const TEN& M)
{
    return M*2;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r+s==1,SYMMETRIC_TENSOR<T,3-r-s,m,n> >::type
Twice_Symmetric_Part(const TENSOR<T,n,n,m>& M)
{
    SYMMETRIC_TENSOR<T,3-r-s,m,n> sm;
    for(int k=0;k<m;k++)
        for(int i=0;i<n;i++)
            for(int j=i;j<n;j++)
                sm.x(k)(i,j)=M.x(i)(j,k)+M.x(j)(i,k);
    return sm;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r+s==2,SYMMETRIC_TENSOR<T,3-r-s,m,n> >::type
Twice_Symmetric_Part(const TENSOR<T,n,m,n>& M)
{
    SYMMETRIC_TENSOR<T,3-r-s,m,n> sm;
    for(int j=0;j<m;j++)
        for(int i=0;i<n;i++)
            for(int k=i;k<n;k++)
                sm.x(j)(i,k)=M.x(i)(j,k)+M.x(k)(j,i);
    return sm;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r+s==3,SYMMETRIC_TENSOR<T,3-r-s,m,n> >::type
Twice_Symmetric_Part(const TENSOR<T,m,n,n>& M)
{
    SYMMETRIC_TENSOR<T,3-r-s,m,n> sm;
    for(int i=0;i<m;i++)
        sm.x(i)=M.x(i).Twice_Symmetric_Part();
    return sm;
}

template<int r,int s,int u,class T,int m>
typename enable_if<3-r-s!=u,SYMMETRIC_TENSOR<T,3-r-s,m,m> >::type
Twice_Symmetric_Part(const SYMMETRIC_TENSOR<T,u,m,m>& M)
{
    SYMMETRIC_TENSOR<T,3-r-s,m,m> sm;
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            for(int k=j;k<m;k++)
                sm.x(i)(j,k)=M.x(j)(i,k)+M.x(k)(i,j);
    return sm;
}

template<int r,int s,int u,class T,int m>
typename enable_if<3-r-s!=u,SYMMETRIC_TENSOR<T,3-r-s,m,m> >::type
Twice_Symmetric_Part(const VEC_ID_SYM_TENSOR<T,u,m>& M)
{
    return VEC_ID_SYM_TENSOR<T,3-r-s,m>(M.v)+VEC_ID_TENSOR<T,3-r-s,m,m>(M.v*2);
}

// template<int r,int s,class T,int m,int n,int p>
// typename enable_if<r+s==1?m==n:r+s==2?m==p:n==p,ZERO_TENSOR<T,m,n,p> >::type
// Twice_Symmetric_Part(const ZERO_TENSOR<T,m,n,p>& M)
// {
//     return M;
// }

template<int r,int s,class T,int m,int n> VEC_ID_TENSOR<T,3-r-s,m,n>
Twice_Symmetric_Part(const VEC_ID_TENSOR<T,3-r-s,m,n>& M)
{
    return M*2;
}

template<int r,int s,class T,int u,int m>
typename enable_if<u!=3-r-s,VEC_ID_SYM_TENSOR<T,3-r-s,m> >::type
Twice_Symmetric_Part(const VEC_ID_TENSOR<T,u,m,m>& M)
{
    return VEC_ID_SYM_TENSOR<T,3-r-s,m>(M.v);
}

template<int r,int s,class T> ZERO_TENSOR<T,3>
Twice_Symmetric_Part(const PERMUTATION_TENSOR<T>& M)
{
    return ZERO_TENSOR<T,3>();
}

template<int r,int s,class MAT,class VEC>
typename enable_if<(IS_MATRIX<MAT>::value && IS_VECTOR<VEC>::value),decltype(Twice_Symmetric_Part<r,s>(Tensor_Product<r>(MAT(),VEC())))>::type
Symmetric_Tensor_Product(const MAT& a,const VEC& b)
{
    return Twice_Symmetric_Part<r,s>(Tensor_Product<r>(a,b));
}

//#####################################################################
// operator+=
//##################################################################### 

template<class T,int m,int n,int p> TENSOR<T,m,n,p>&
operator+=(TENSOR<T,m,n,p>& a,const TENSOR<T,m,n,p>& b)
{
    a.x+=b.x;
    return a;
}

template<class T,int m,int n,int p> TENSOR<T,m,n,p>&
operator+=(TENSOR<T,m,n,p>& a,const ZERO_TENSOR<T,m,n,p>& b)
{
    return a;
}

template<class T,int m> TENSOR<T,m>&
operator+=(TENSOR<T,m>& a,const DIAGONAL_TENSOR<T,m>& b)
{
    for(int i=0;i<m;i++) a.x(i)(i,i)+=b.v(i);
    return a;
}

template<class T,int m,int n> TENSOR<T,m,n,n>&
operator+=(TENSOR<T,m,n,n>& a,const SYMMETRIC_TENSOR<T,0,m,n>& b)
{
    for(int i=0;i<m;i++) a.x(i)+=b.x(i);
    return a;
}

template<class T,int m,int n> TENSOR<T,n,m,n>&
operator+=(TENSOR<T,n,m,n>& a,const SYMMETRIC_TENSOR<T,1,m,n>& b)
{
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) a.x(j).Add_Row(i,b.x(i).Row(j));
    return a;
}

template<class T,int m,int n> TENSOR<T,n,n,m>&
operator+=(TENSOR<T,n,n,m>& a,const SYMMETRIC_TENSOR<T,2,m,n>& b)
{
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) a.x(j).Add_Column(i,b.x(i).Column(j));
    return a;
}

template<class T,int m,int n> TENSOR<T,m,n,n>&
operator+=(TENSOR<T,m,n,n>& a,const VEC_ID_TENSOR<T,0,m,n>& b)
{
    for(int i=0;i<m;i++) a.x(i)+=b.v(i);
    return a;
}

template<class T,int m,int n> TENSOR<T,n,m,n>&
operator+=(TENSOR<T,n,m,n>& a,const VEC_ID_TENSOR<T,1,m,n>& b)
{
    for(int i=0;i<n;i++) a.x(i).Add_Column(i,b.v);
    return a;
}

template<class T,int m,int n> TENSOR<T,n,n,m>&
operator+=(TENSOR<T,n,n,m>& a,const VEC_ID_TENSOR<T,2,m,n>& b)
{
    for(int i=0;i<n;i++) a.x(i).Add_Row(i,b.v);
    return a;
}

template<class T,int u,int m> TENSOR<T,m>
operator+=(TENSOR<T,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    a+=VEC_ID_TENSOR<T,u==0?1:0,m>(b.v);
    a+=VEC_ID_TENSOR<T,u==2?1:2,m>(b.v);
    return a;
}

template<class T> TENSOR<T,3,3,3>&
operator+=(TENSOR<T,3,3,3>& a,const PERMUTATION_TENSOR<T>& b)
{
    a.x(0)(1,2)+=b.x;
    a.x(1)(2,0)+=b.x;
    a.x(2)(0,1)+=b.x;
    a.x(0)(2,1)-=b.x;
    a.x(1)(0,2)-=b.x;
    a.x(2)(1,0)-=b.x;
    return a;
}

template<class T,int u,int m,int n> SYMMETRIC_TENSOR<T,u,m,n>&
operator+=(SYMMETRIC_TENSOR<T,u,m,n>& a,const SYMMETRIC_TENSOR<T,u,m,n>& b)
{
    a.x+=b.x;
    return a;
}

template<class T,int u,int m,int n> SYMMETRIC_TENSOR<T,u,m,n>&
operator+=(SYMMETRIC_TENSOR<T,u,m,n>& a,const ZERO_TENSOR<T,u==0?m:n,u==1?m:n,u==2?m:n>& b)
{
    return a;
}

template<class T,int u,int m,int n> SYMMETRIC_TENSOR<T,u,m,n>&
operator+=(SYMMETRIC_TENSOR<T,u,m,n>& a,const VEC_ID_TENSOR<T,u,m,n>& b)
{
    for(int i=0;i<m;i++) a.x(i)+=b.v(i);
    return a;
}

template<class T,int u,int m> SYMMETRIC_TENSOR<T,u,m>
operator+=(SYMMETRIC_TENSOR<T,u,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    TENSOR<T,m> t;
    t+=VEC_ID_TENSOR<T,1,m,m>(b.v);
    for(int i=0;i<m;i++) a.x(i)+=t.x(i).Twice_Symmetric_Part();
    return a;
}

template<class T,int u,int m> SYMMETRIC_TENSOR<T,u,m,m>&
operator+=(SYMMETRIC_TENSOR<T,u,m,m>& a,const DIAGONAL_TENSOR<T,m>& b)
{
    for(int i=0;i<m;i++) a.x(i)(i,i)+=b.v(i);
    return a;
}

//#####################################################################
// operator-=
//##################################################################### 

template<class T,int m,int n,int p> TENSOR<T,m,n,p>&
operator-=(TENSOR<T,m,n,p>& a,const TENSOR<T,m,n,p>& b)
{
    a.x-=b.x;
    return a;
}

template<class T,int m,int n,int p> TENSOR<T,m,n,p>&
operator-=(TENSOR<T,m,n,p>& a,const ZERO_TENSOR<T,m,n,p>& b)
{
    return a;
}

template<class T,int m,int n> TENSOR<T,m,n,n>&
operator-=(TENSOR<T,m,n,n>& a,const SYMMETRIC_TENSOR<T,0,m,n>& b)
{
    for(int i=0;i<m;i++) a.x(i)-=b.x(i);
    return a;
}

template<class T,int m,int n> TENSOR<T,n,m,n>&
operator-=(TENSOR<T,n,m,n>& a,const SYMMETRIC_TENSOR<T,1,m,n>& b)
{
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) a.x(j).Add_Row(i,-b.x(i).Row(j));
    return a;
}

template<class T,int m,int n> TENSOR<T,n,n,m>&
operator-=(TENSOR<T,n,n,m>& a,const SYMMETRIC_TENSOR<T,2,m,n>& b)
{
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) a.x(j).Add_Column(i,-b.x(i).Column(j));
    return a;
}

template<class T,int m,int n> TENSOR<T,m,n,n>&
operator-=(TENSOR<T,m,n,n>& a,const VEC_ID_TENSOR<T,0,m,n>& b)
{
    for(int i=0;i<m;i++) a.x(i)-=b.v(i);
    return a;
}

template<class T,int m,int n> TENSOR<T,n,m,n>&
operator-=(TENSOR<T,n,m,n>& a,const VEC_ID_TENSOR<T,1,m,n>& b)
{
    return a+=-b;
}

template<class T,int m,int n> TENSOR<T,n,n,m>&
operator-=(TENSOR<T,n,n,m>& a,const VEC_ID_TENSOR<T,2,m,n>& b)
{
    return a+=-b;
}

template<class T,int u,int m> TENSOR<T,m>
operator-=(TENSOR<T,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    return a+=-b;
}

template<class T> TENSOR<T,3>&
operator-=(TENSOR<T,3>& a,const PERMUTATION_TENSOR<T>& b)
{
    a.x(0)(1,2)-=b.x;
    a.x(1)(2,0)-=b.x;
    a.x(2)(0,1)-=b.x;
    a.x(0)(2,1)+=b.x;
    a.x(1)(0,2)+=b.x;
    a.x(2)(1,0)+=b.x;
    return a;
}

template<class T,int m> TENSOR<T,m>&
operator-=(TENSOR<T,m>& a,const DIAGONAL_TENSOR<T,m>& b)
{
    for(int i=0;i<m;i++) a.x(i)(i,i)-=b.v(i);
    return a;
}

template<class T,int u,int m,int n> SYMMETRIC_TENSOR<T,u,m,n>&
operator-=(SYMMETRIC_TENSOR<T,u,m,n>& a,const SYMMETRIC_TENSOR<T,u,m,n>& b)
{
    a.x-=b.x;
    return a;
}

template<class T,int u,int m,int n> SYMMETRIC_TENSOR<T,u,m,n>&
operator-=(SYMMETRIC_TENSOR<T,u,m,n>& a,const ZERO_TENSOR<T,u==0?m:n,u==1?m:n,u==2?m:n>& b)
{
    return a;
}

template<class T,int u,int m,int n> SYMMETRIC_TENSOR<T,u,m,n>&
operator-=(SYMMETRIC_TENSOR<T,u,m,n>& a,const VEC_ID_TENSOR<T,u,m,n>& b)
{
    for(int i=0;i<m;i++) a.x(i)-=b.v(i);
    return a;
}

template<class T,int u,int m> SYMMETRIC_TENSOR<T,u,m>
operator-=(SYMMETRIC_TENSOR<T,u,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    TENSOR<T,m> t;
    t+=VEC_ID_TENSOR<T,1,m,m>(b.v);
    for(int i=0;i<m;i++) a.x(i)-=t.x(i).Twice_Symmetric_Part();
    return a;
}

template<class T,int u,int m> SYMMETRIC_TENSOR<T,u,m,m>&
operator-=(SYMMETRIC_TENSOR<T,u,m,m>& a,const DIAGONAL_TENSOR<T,m>& b)
{
    for(int i=0;i<m;i++) a.x(i)(i,i)-=b.v(i);
    return a;
}

//#####################################################################
// operator+
//##################################################################### 

template<class T> struct TENSOR_ORDER {};
template<class T,int m,int n,int p> struct TENSOR_ORDER<ZERO_TENSOR<T,m,n,p> > {static const int value=0;};
template<class T,int m,int n,int p> struct TENSOR_ORDER<TENSOR<T,m,n,p> > {static const int value=1;};
template<class T,int m,int n,int u> struct TENSOR_ORDER<SYMMETRIC_TENSOR<T,u,m,n> > {static const int value=2;};
template<class T,int m,int n,int u> struct TENSOR_ORDER<VEC_ID_TENSOR<T,u,m,n> > {static const int value=3;};
template<class T,int u,int m> struct TENSOR_ORDER<VEC_ID_SYM_TENSOR<T,u,m> > {static const int value=4;};
template<class T> struct TENSOR_ORDER<PERMUTATION_TENSOR<T> > {static const int value=5;};
template<class T,int m> struct TENSOR_ORDER<DIAGONAL_TENSOR<T,m> > {static const int value=6;};

template<class TEN,class T> TEN
operator+(const ZERO_TENSOR<T,TEN::m,TEN::n,TEN::p>& a,const TEN& b)
{
    return b;
}

template<class TEN,class T>
typename enable_if<(0<TENSOR_ORDER<TEN>::value),TEN>::type
operator+(const TEN& b,const ZERO_TENSOR<T,TEN::m,TEN::n,TEN::p>& a)
{
    return b;
}

template<class TEN,class T>
typename enable_if<(1<=TENSOR_ORDER<TEN>::value),TENSOR<T,TEN::m,TEN::n,TEN::p> >::type
operator+(const TENSOR<T,TEN::m,TEN::n,TEN::p>& a,const TEN& b)
{
    TENSOR<T,TEN::m,TEN::n,TEN::p> t(a);
    return t+=b;
}

template<class TEN,class T>
typename enable_if<(1<TENSOR_ORDER<TEN>::value),TENSOR<T,TEN::m,TEN::n,TEN::p> >::type
operator+(const TEN& b,const TENSOR<T,TEN::m,TEN::n,TEN::p>& a)
{
    TENSOR<T,TEN::m,TEN::n,TEN::p> t(a);
    return t+=b;
}

template<class TEN,class T,int u,int m,int n>
typename enable_if<(2<=TENSOR_ORDER<TEN>::value && IS_SYM_TENSOR<u,TEN>::value),SYMMETRIC_TENSOR<T,u,m,n> >::type
operator+(const SYMMETRIC_TENSOR<T,u,m,n>& a,const TEN& b)
{
    SYMMETRIC_TENSOR<T,u,m,n> t(a);
    return t+=b;
}

template<class TEN,class T,int u,int m,int n>
typename enable_if<(2<TENSOR_ORDER<TEN>::value && IS_SYM_TENSOR<u,TEN>::value),SYMMETRIC_TENSOR<T,u,m,n> >::type
operator+(const TEN& b,const SYMMETRIC_TENSOR<T,u,m,n>& a)
{
    SYMMETRIC_TENSOR<T,u,m,n> t(a);
    return t+=b;
}

template<class TEN,class T,int u,int m,int n>
typename enable_if<(2<=TENSOR_ORDER<TEN>::value && !IS_SYM_TENSOR<u,TEN>::value),TENSOR<T,TEN::m,TEN::n,TEN::p> >::type
operator+(const SYMMETRIC_TENSOR<T,u,m,n>& a,const TEN& b)
{
    TENSOR<T,TEN::m,TEN::n,TEN::p> t;
    t+=a;
    return t+=b;
}

template<class TEN,class T,int u,int m,int n>
typename enable_if<(2<TENSOR_ORDER<TEN>::value && !IS_SYM_TENSOR<u,TEN>::value),TENSOR<T,u==0?m:n,u==1?m:n,u==2?m:n> >::type
operator+(const TEN& b,const SYMMETRIC_TENSOR<T,u,m,n>& a)
{
    TENSOR<T,TEN::m,TEN::n,TEN::p> t;
    t+=a;
    return t+=b;
}

template<class T,int u,int m,int n> VEC_ID_TENSOR<T,u,m,n>
operator+(const VEC_ID_TENSOR<T,u,m,n>& a,const VEC_ID_TENSOR<T,u,m,n>& b)
{
    return VEC_ID_TENSOR<T,u,m,n>(a.v+b.v);
}

template<class T,int u,int m,int n> SYMMETRIC_TENSOR<T,u,m>
operator+(const VEC_ID_TENSOR<T,u,m,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    SYMMETRIC_TENSOR<T,u,m> t;
    t+=a;
    t+=b;
    return t;
}

template<class T,int u,int m> SYMMETRIC_TENSOR<T,u,m>
operator+(const VEC_ID_SYM_TENSOR<T,u,m>& b,const VEC_ID_TENSOR<T,u,m,m>& a)
{
    SYMMETRIC_TENSOR<T,u,m> t;
    t+=a;
    t+=b;
    return t;
}

template<class T,int u,int m> SYMMETRIC_TENSOR<T,u,m>
operator+(const VEC_ID_TENSOR<T,u,m,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    SYMMETRIC_TENSOR<T,u,m> t;
    t+=a;
    t+=b;
    return t;
}

template<class T,int u,int m> VEC_ID_SYM_TENSOR<T,u,m>
operator+(const VEC_ID_SYM_TENSOR<T,u,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    return VEC_ID_SYM_TENSOR<T,u,m>(a.v+b.v);
}

template<class T,int u,int v,int m> TENSOR<T,m>
operator+(const VEC_ID_SYM_TENSOR<T,u,m>& a,const VEC_ID_SYM_TENSOR<T,v,m>& b)
{
    TENSOR<T,m> t;
    t+=a;
    return t+=b;
}

template<class T,int u,int v,int m> TENSOR<T,m>
operator+(const VEC_ID_SYM_TENSOR<T,u,m>& a,const VEC_ID_TENSOR<T,v,m,m>& b)
{
    TENSOR<T,m> t;
    t+=a;
    return t+=b;
}

template<class T,int u,int v,int m> TENSOR<T,m>
operator+(const VEC_ID_TENSOR<T,v,m,m>& b,const VEC_ID_SYM_TENSOR<T,u,m>& a)
{
    TENSOR<T,m> t;
    t+=a;
    return t+=b;
}

template<class T,int u,int v,int m> typename enable_if<u!=v,TENSOR<T,m> >::type
operator+(const VEC_ID_TENSOR<T,v,m,m>& b,const VEC_ID_TENSOR<T,u,m,m>& a)
{
    TENSOR<T,m> t;
    t+=a;
    return t+=b;
}

template<class T> PERMUTATION_TENSOR<T>
operator+(const PERMUTATION_TENSOR<T>& a,const PERMUTATION_TENSOR<T>& b)
{
    return PERMUTATION_TENSOR<T>(a.x+b.x);
}

template<class T,int u> TENSOR<T,3,3,3>
operator+(const PERMUTATION_TENSOR<T>& a,const SYMMETRIC_TENSOR<T,u,3,3>& b)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t+=b;
    return t;
}

template<class T,int u> TENSOR<T,3,3,3>
operator+(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR<T,u,3,3>& b)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t+=b;
    return t;
}

template<class T,int u> TENSOR<T,3,3,3>
operator+(const PERMUTATION_TENSOR<T>& a,const VEC_ID_SYM_TENSOR<T,u,3>& b)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t+=b;
    return t;
}

template<class T> TENSOR<T,3,3,3>
operator+(const PERMUTATION_TENSOR<T>& a,const TENSOR<T,3,3,3>& b)
{
    TENSOR<T,3,3,3> t=b;
    return t+=a;
}

template<class T,int u> TENSOR<T,3,3,3>
operator+(const SYMMETRIC_TENSOR<T,u,3,3>& b,const PERMUTATION_TENSOR<T>& a)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t+=b;
    return t;
}

template<class T,int u> TENSOR<T,3,3,3>
operator+(const VEC_ID_TENSOR<T,u,3,3>& b,const PERMUTATION_TENSOR<T>& a)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t+=b;
    return t;
}

template<class T,int u> TENSOR<T,3,3,3>
operator+(const VEC_ID_SYM_TENSOR<T,u,3>& b,const PERMUTATION_TENSOR<T>& a)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t+=b;
    return t;
}

template<class T> TENSOR<T,3,3,3>
operator+(const TENSOR<T,3,3,3>& b,const PERMUTATION_TENSOR<T>& a)
{
    TENSOR<T,3,3,3> t=b;
    return t+=a;
}

//#####################################################################
// operator-
//##################################################################### 

template<class TEN,class T> TEN
operator-(const TEN& b,const ZERO_TENSOR<T,TEN::m,TEN::n,TEN::p>& a)
{
    return b;
}

template<class TEN,class T>
typename enable_if<TENSOR_ORDER<TEN>::value!=0,TEN>::type
operator-(const ZERO_TENSOR<T,TEN::m,TEN::n,TEN::p>& a,const TEN& b)
{
    return -b;
}

template<class TEN,class T>
typename enable_if<(1<=TENSOR_ORDER<TEN>::value),TENSOR<T,TEN::m,TEN::n,TEN::p> >::type
operator-(const TENSOR<T,TEN::m,TEN::n,TEN::p>& a,const TEN& b)
{
    TENSOR<T,TEN::m,TEN::n,TEN::p> t(a);
    return t-=b;
}

template<class TEN,class T>
typename enable_if<(2<=TENSOR_ORDER<TEN>::value),TENSOR<T,TEN::m,TEN::n,TEN::p> >::type
operator-(const TEN& b,const TENSOR<T,TEN::m,TEN::n,TEN::p>& a)
{
    TENSOR<T,TEN::m,TEN::n,TEN::p> t(-a);
    return t+=b;
}

template<class TEN,class T,int u,int m,int n>
typename enable_if<(2<=TENSOR_ORDER<TEN>::value && IS_SYM_TENSOR<u,TEN>::value),SYMMETRIC_TENSOR<T,u,m,n> >::type
operator-(const SYMMETRIC_TENSOR<T,u,m,n>& a,const TEN& b)
{
    SYMMETRIC_TENSOR<T,u,m,n> t(a);
    return t-=b;
}

template<class TEN,class T,int u,int m,int n>
typename enable_if<(2<TENSOR_ORDER<TEN>::value && IS_SYM_TENSOR<u,TEN>::value),SYMMETRIC_TENSOR<T,u,m,n> >::type
operator-(const TEN& b,const SYMMETRIC_TENSOR<T,u,m,n>& a)
{
    SYMMETRIC_TENSOR<T,u,m,n> t(-a);
    return t+=b;
}

template<class TEN,class T,int u,int m,int n>
typename enable_if<(2<=TENSOR_ORDER<TEN>::value && !IS_SYM_TENSOR<u,TEN>::value),TENSOR<T,u==0?m:n,u==1?m:n,u==2?m:n> >::type
operator-(const SYMMETRIC_TENSOR<T,u,m,n>& a,const TEN& b)
{
    TENSOR<T,u==0?m:n,u==1?m:n,u==2?m:n> t;
    t+=a;
    return t-=b;
}

template<class TEN,class T,int u,int m,int n>
typename enable_if<(2<TENSOR_ORDER<TEN>::value && !IS_SYM_TENSOR<u,TEN>::value),TENSOR<T,u==0?m:n,u==1?m:n,u==2?m:n> >::type
operator-(const TEN& b,const SYMMETRIC_TENSOR<T,u,m,n>& a)
{
    TENSOR<T,u==0?m:n,u==1?m:n,u==2?m:n> t;
    t-=a;
    return t+=b;
}

template<class T,int u,int m,int n> VEC_ID_TENSOR<T,u,m,n>
operator-(const VEC_ID_TENSOR<T,u,m,n>& a,const VEC_ID_TENSOR<T,u,m,n>& b)
{
    return VEC_ID_TENSOR<T,u,m,n>(a.v-b.v);
}

template<class T,int u,int m,int n> SYMMETRIC_TENSOR<T,u,m>
operator-(const VEC_ID_TENSOR<T,u,m,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    SYMMETRIC_TENSOR<T,u,m> t;
    t+=a;
    t-=b;
    return t;
}

template<class T,int u,int m> SYMMETRIC_TENSOR<T,u,m>
operator-(const VEC_ID_SYM_TENSOR<T,u,m>& a,const VEC_ID_TENSOR<T,u,m,m>& b)
{
    SYMMETRIC_TENSOR<T,u,m> t;
    t+=a;
    t-=b;
    return t;
}

template<class T,int u,int v,int m> TENSOR<T,m>
operator-(const VEC_ID_SYM_TENSOR<T,u,m>& a,const VEC_ID_TENSOR<T,v,m,m>& b)
{
    TENSOR<T,m> t;
    t+=a;
    return t-=b;
}

template<class T,int u,int v,int m> TENSOR<T,m>
operator-(const VEC_ID_TENSOR<T,v,m,m>& b,const VEC_ID_SYM_TENSOR<T,u,m>& a)
{
    TENSOR<T,m> t;
    t-=a;
    return t+=b;
}

template<class T,int u,int m> SYMMETRIC_TENSOR<T,u,m>
operator-(const VEC_ID_TENSOR<T,u,m,m>& b,const VEC_ID_SYM_TENSOR<T,u,m>& a)
{
    SYMMETRIC_TENSOR<T,u,m> t;
    t-=a;
    t+=b;
    return t;
}

template<class T,int u,int v,int m> typename enable_if<u!=v,TENSOR<T,m> >::type
operator-(const VEC_ID_TENSOR<T,v,m,m>& b,const VEC_ID_TENSOR<T,u,m,m>& a)
{
    TENSOR<T,m> t;
    t-=a;
    return t+=b;
}

template<class T,int u,int m> VEC_ID_SYM_TENSOR<T,u,m>
operator-(const VEC_ID_SYM_TENSOR<T,u,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    return VEC_ID_SYM_TENSOR<T,u,m>(a.v-b.v);
}

template<class TEN0,class TEN1>
typename enable_if<(TENSOR_ORDER<TEN0>::value<=TENSOR_ORDER<TEN1>::value),TENSOR<typename TEN0::SCALAR,TEN0::m,TEN0::n,TEN0::p> >::type
operator-(const TEN0& a,const TEN1& b)
{
    TENSOR<typename TEN0::SCALAR,TEN0::m,TEN0::n,TEN0::p> t;
    t+=a;
    t-=b;
    return t;
}

template<class T> PERMUTATION_TENSOR<T>
operator-(const PERMUTATION_TENSOR<T>& a,const PERMUTATION_TENSOR<T>& b)
{
    return PERMUTATION_TENSOR<T>(a.x-b.x);
}

template<class T,int u> TENSOR<T,3,3,3>
operator-(const PERMUTATION_TENSOR<T>& a,const SYMMETRIC_TENSOR<T,u,3,3>& b)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t-=b;
    return t;
}

template<class T,int u> TENSOR<T,3,3,3>
operator-(const PERMUTATION_TENSOR<T>& a,const VEC_ID_TENSOR<T,u,3,3>& b)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t-=b;
    return t;
}

template<class T,int u> TENSOR<T,3,3,3>
operator-(const PERMUTATION_TENSOR<T>& a,const VEC_ID_SYM_TENSOR<T,u,3>& b)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t-=b;
    return t;
}

template<class T> TENSOR<T,3,3,3>
operator-(const PERMUTATION_TENSOR<T>& a,const TENSOR<T,3,3,3>& b)
{
    TENSOR<T,3,3,3> t;
    t+=a;
    t-=b;
    return t;
}

//#####################################################################
// Function Contract<r,s>(T,M)
//##################################################################### 

template<int r,int s,class TEN,class T,int q0,int q1>
typename enable_if<IS_TENSOR<TEN>::value && (r==0?TEN::m:r==1?TEN::n:TEN::p)==(s==0?q0:q1),ZERO_TENSOR<T,r==0?TEN::m:(s==0?q0:q1),r==1?TEN::n:(s==0?q0:q1),r==2?TEN::p:(s==0?q0:q1)> >::type
Contract(const TEN& a,const ZERO_MATRIX<T,q0,q1>& M)
{
    return ZERO_TENSOR<T,r==0?TEN::m:(s==0?q0:q1),r==1?TEN::n:(s==0?q0:q1),r==2?TEN::p:(s==0?q0:q1)>();
}

template<int r,int s,class TEN,class T>
typename enable_if<IS_TENSOR<TEN>::value,TEN>::type
Contract(const TEN& a,const SCALE_MATRIX<T,r==0?TEN::m:r==1?TEN::n:TEN::p>& M)
{
    return a*M.x;
}

template<int r,int s,class TEN,class T>
typename enable_if<IS_TENSOR<TEN>::value,TEN>::type
Contract(const TEN& a,const IDENTITY_MATRIX<T,r==0?TEN::m:r==1?TEN::n:TEN::p>& M)
{
    return a;
}

template<int r,int s,class T,int m,int n,int p,int q0,int q1>
typename enable_if<((r==0?m:r==1?n:p)==(s==0?q0:q1)),ZERO_TENSOR<T,r!=0?m:(s==1?q0:q1),r!=1?n:(s==1?q0:q1),r!=2?p:(s==1?q0:q1)> >::type
Contract(const ZERO_TENSOR<T,m,n,p>& a,const MATRIX<T,q0,q1>& M)
{
    return ZERO_TENSOR<T,r!=0?m:(s==1?q0:q1),r!=1?n:(s==1?q0:q1),r!=2?p:(s==1?q0:q1)>();
}

template<int r,int s,class T,int m,int n,int p> ZERO_TENSOR<T,m,n,p>
Contract(const ZERO_TENSOR<T,m,n,p>& a,const SYMMETRIC_MATRIX<T,r==0?m:r==1?n:p>& M)
{
    return ZERO_TENSOR<T,m,n,p>();
}

template<int r,int s,class T,int m,int n,int p> ZERO_TENSOR<T,m,n,p>
Contract(const ZERO_TENSOR<T,m,n,p>& a,const DIAGONAL_MATRIX<T,r==0?m:r==1?n:p>& M)
{
    return ZERO_TENSOR<T,m,n,p>();
}

template<int r,int s,class T,int m,int n,int p,int q>
typename enable_if<r==0 && s==0,TENSOR<T,q,n,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const MATRIX<T,m,q>& M)
{
    TENSOR<T,q,n,p> t;
    for(int i=0;i<m;i++)
        for(int l=0;l<q;l++)
            t.x(l)+=a.x(i)*M(i,l);
    return t;
}

template<int r,int s,class T,int m,int n,int p,int q>
typename enable_if<r==1 && s==0,TENSOR<T,m,q,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const MATRIX<T,n,q>& M)
{
    TENSOR<T,m,q,p> t;
    for(int i=0;i<m;i++)
        t.x(i)+=M.Transpose_Times(a.x(i));
    return t;
}

template<int r,int s,class T,int m,int n,int p,int q>
typename enable_if<r==2 && s==0,TENSOR<T,m,n,q> >::type
Contract(const TENSOR<T,m,n,p>& a,const MATRIX<T,p,q>& M)
{
    TENSOR<T,m,n,q> t;
    for(int i=0;i<m;i++)
        t.x(i)+=a.x(i)*M;
    return t;
}

template<int r,int s,class T,int m,int n,int p,int q>
typename enable_if<r==0 && s==1,TENSOR<T,q,n,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const MATRIX<T,q,m>& M)
{
    TENSOR<T,q,n,p> t;
    for(int i=0;i<m;i++)
        for(int l=0;l<q;l++)
            t.x(l)+=a.x(i)*M(l,i);
    return t;
}

template<int r,int s,class T,int m,int n,int p,int q>
typename enable_if<r==1 && s==1,TENSOR<T,m,q,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const MATRIX<T,q,n>& M)
{
    TENSOR<T,m,q,p> t;
    for(int i=0;i<m;i++)
        t.x(i)+=M*a.x(i);
    return t;
}

template<int r,int s,class T,int m,int n,int p,int q>
typename enable_if<r==2 && s==1,TENSOR<T,m,n,q> >::type
Contract(const TENSOR<T,m,n,p>& a,const MATRIX<T,q,p>& M)
{
    TENSOR<T,m,n,q> t;
    for(int i=0;i<m;i++)
        t.x(i)+=a.x(i).Times_Transpose(M);
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<r==0,TENSOR<T,m,n,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const SYMMETRIC_MATRIX<T,m>& M)
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<m;i++)
        for(int l=0;l<m;l++)
            t.x(l)+=a.x(i)*M(i,l);
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<r==1,TENSOR<T,m,n,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const SYMMETRIC_MATRIX<T,n>& M)
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<m;i++)
        t.x(i)+=M.Transpose_Times(a.x(i));
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<r==2,TENSOR<T,m,n,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const SYMMETRIC_MATRIX<T,p>& M)
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<m;i++)
        t.x(i)+=a.x(i)*M;
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<r==0,TENSOR<T,m,n,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const DIAGONAL_MATRIX<T,m>& M)
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<m;i++)
        t.x(i)+=a.x(i)*M(i,i);
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<r==1,TENSOR<T,m,n,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const DIAGONAL_MATRIX<T,n>& M)
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<m;i++)
        t.x(i)+=M*a.x(i);
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<r==2,TENSOR<T,m,n,p> >::type
Contract(const TENSOR<T,m,n,p>& a,const DIAGONAL_MATRIX<T,p>& M)
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<m;i++)
        t.x(i)+=a.x(i)*M;
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<s==0,SYMMETRIC_TENSOR<T,r,p,n> >::type
Contract(const SYMMETRIC_TENSOR<T,r,m,n>& a,const MATRIX<T,m,p>& M)
{
    SYMMETRIC_TENSOR<T,r,p,n> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<p;j++)
            t.x(j)+=a.x(i)*M(i,j);
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==1 && s==0),TENSOR<T,m,p,n> >::type
Contract(const SYMMETRIC_TENSOR<T,0,m,n>& a,const MATRIX<T,n,p>& M)
{
    TENSOR<T,m,p,n> t;
    for(int i=0;i<m;i++)
        t.x(i)=M.Transpose_Times(a.x(i));
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==2 && s==0),TENSOR<T,m,n,p> >::type
Contract(const SYMMETRIC_TENSOR<T,0,m,n>& a,const MATRIX<T,n,p>& M)
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<m;i++)
        t.x(i)=a.x(i)*M;
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==0 && s==0),TENSOR<T,p,m,n> >::type
Contract(const SYMMETRIC_TENSOR<T,1,m,n>& a,const MATRIX<T,n,p>& M)
{
    TENSOR<T,p,m,n> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<p;j++)
            t.x(j).Set_Row(i,a.x(i).Transpose_Times(M.Column(j)));
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==2 && s==0),TENSOR<T,n,m,p> >::type
Contract(const SYMMETRIC_TENSOR<T,1,m,n>& a,const MATRIX<T,n,p>& M)
{
    TENSOR<T,n,m,p> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            t.x(j).Set_Row(i,M.Transpose_Times(a.x(i).Row(j)));
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==0 && s==0),TENSOR<T,p,n,m> >::type
Contract(const SYMMETRIC_TENSOR<T,2,m,n>& a,const MATRIX<T,n,p>& M)
{
    TENSOR<T,p,n,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<p;j++)
            t.x(j).Set_Column(i,a.x(i).Transpose_Times(M.Column(j)));
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==1 && s==0),TENSOR<T,n,p,m> >::type
Contract(const SYMMETRIC_TENSOR<T,2,m,n>& a,const MATRIX<T,n,p>& M)
{
    TENSOR<T,n,p,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            t.x(j).Set_Column(i,M.Transpose_Times(a.x(i).Row(j)));
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<s==1,SYMMETRIC_TENSOR<T,r,p,n> >::type
Contract(const SYMMETRIC_TENSOR<T,r,m,n>& a,const MATRIX<T,p,m>& M)
{
    SYMMETRIC_TENSOR<T,r,p,n> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<p;j++)
            t.x(j)+=a.x(i)*M(j,i);
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==1 && s==1),TENSOR<T,m,p,n> >::type
Contract(const SYMMETRIC_TENSOR<T,0,m,n>& a,const MATRIX<T,p,n>& M)
{
    TENSOR<T,m,p,n> t;
    for(int i=0;i<m;i++)
        t.x(i)=M*a.x(i);
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==2 && s==1),TENSOR<T,m,n,p> >::type
Contract(const SYMMETRIC_TENSOR<T,0,m,n>& a,const MATRIX<T,p,n>& M)
{
    TENSOR<T,m,n,p> t;
    for(int i=0;i<m;i++)
        t.x(i)=a.x(i).Times_Transpose(M);
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==0 && s==1),TENSOR<T,p,m,n> >::type
Contract(const SYMMETRIC_TENSOR<T,1,m,n>& a,const MATRIX<T,p,n>& M)
{
    TENSOR<T,p,m,n> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<p;j++)
            t.x(j).Set_Row(i,a.x(i).Transpose_Times(M.Row(j)));
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==2 && s==1),TENSOR<T,n,m,p> >::type
Contract(const SYMMETRIC_TENSOR<T,1,m,n>& a,const MATRIX<T,p,n>& M)
{
    TENSOR<T,n,m,p> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            t.x(j).Set_Row(i,M*a.x(i).Row(j));
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==0 && s==1),TENSOR<T,p,n,m> >::type
Contract(const SYMMETRIC_TENSOR<T,2,m,n>& a,const MATRIX<T,p,n>& M)
{
    TENSOR<T,p,n,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<p;j++)
            t.x(j).Set_Column(i,a.x(i).Transpose_Times(M.Row(j)));
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<(r==1 && s==1),TENSOR<T,n,p,m> >::type
Contract(const SYMMETRIC_TENSOR<T,2,m,n>& a,const MATRIX<T,p,n>& M)
{
    TENSOR<T,n,p,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            t.x(j).Set_Column(i,M*a.x(i).Row(j));
    return t;
}

template<int r,int s,class T,int m,int n> SYMMETRIC_TENSOR<T,r,m,n>
Contract(const SYMMETRIC_TENSOR<T,r,m,n>& a,const SYMMETRIC_MATRIX<T,m>& M)
{
    SYMMETRIC_TENSOR<T,r,m,n> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            t.x(j)+=a.x(i)*M(i,j);
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==1,TENSOR<T,m,n,n> >::type
Contract(const SYMMETRIC_TENSOR<T,0,m,n>& a,const SYMMETRIC_MATRIX<T,n>& M)
{
    TENSOR<T,m,n,n> t;
    for(int i=0;i<m;i++)
        t.x(i)=M.Transpose_Times(a.x(i));
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==2,TENSOR<T,m,n,n> >::type
Contract(const SYMMETRIC_TENSOR<T,0,m,n>& a,const SYMMETRIC_MATRIX<T,n>& M)
{
    TENSOR<T,m,n,n> t;
    for(int i=0;i<m;i++)
        t.x(i)=a.x(i).Times_Transpose(M);
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==0,TENSOR<T,n,m,n> >::type
Contract(const SYMMETRIC_TENSOR<T,1,m,n>& a,const SYMMETRIC_MATRIX<T,n>& M)
{
    TENSOR<T,n,m,n> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            t.x(j).Set_Row(i,a.x(i).Transpose_Times(M.Column(j)));
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==2,TENSOR<T,n,m,n> >::type
Contract(const SYMMETRIC_TENSOR<T,1,m,n>& a,const SYMMETRIC_MATRIX<T,n>& M)
{
    TENSOR<T,n,m,n> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            t.x(j).Set_Row(i,M*a.x(i).Row(j));
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==0,TENSOR<T,n,n,m> >::type
Contract(const SYMMETRIC_TENSOR<T,2,m,n>& a,const SYMMETRIC_MATRIX<T,n>& M)
{
    TENSOR<T,n,n,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            t.x(j).Set_Column(i,a.x(i).Transpose_Times(M.Column(j)));
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==1,TENSOR<T,n,n,m> >::type
Contract(const SYMMETRIC_TENSOR<T,2,m,n>& a,const SYMMETRIC_MATRIX<T,n>& M)
{
    TENSOR<T,n,n,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            t.x(j).Set_Column(i,M.Transpose_Times(a.x(i).Row(j)));
    return t;
}

template<int r,int s,class T,int m> SYMMETRIC_TENSOR<T,r,m>
Contract(const DIAGONAL_TENSOR<T,m>& a,const SYMMETRIC_MATRIX<T,m>& M)
{
    SYMMETRIC_TENSOR<T,r,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<m;j++)
            t.x(j)(i,i)+=a.x(i)*M(i,j);
    return t;
}

template<int r,int s,class T,int m,int n> SYMMETRIC_TENSOR<T,r,m,n>
Contract(const SYMMETRIC_TENSOR<T,r,m,n>& a,const DIAGONAL_MATRIX<T,m>& M)
{
    SYMMETRIC_TENSOR<T,r,m,n> t;
    for(int i=0;i<m;i++)
        t.x(i)+=a.x(i)*M(i,i);
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==1,TENSOR<T,m,n,n> >::type
Contract(const SYMMETRIC_TENSOR<T,0,m,n>& a,const DIAGONAL_MATRIX<T,n>& M)
{
    TENSOR<T,m,n,n> t;
    for(int i=0;i<m;i++)
        t.x(i)=M.Transpose_Times(a.x(i));
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==2,TENSOR<T,m,n,n> >::type
Contract(const SYMMETRIC_TENSOR<T,0,m,n>& a,const DIAGONAL_MATRIX<T,n>& M)
{
    TENSOR<T,m,n,n> t;
    for(int i=0;i<m;i++)
        t.x(i)=a.x(i).Times_Transpose(M);
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==0,TENSOR<T,n,m,n> >::type
Contract(const SYMMETRIC_TENSOR<T,1,m,n>& a,const DIAGONAL_MATRIX<T,n>& M)
{
    TENSOR<T,n,m,n> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            for(int k=0;k<n;k++)
                t.x(j)(i,k)+=a.x(i)(j,k)*M(j,j);
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==2,TENSOR<T,n,m,n> >::type
Contract(const SYMMETRIC_TENSOR<T,1,m,n>& a,const DIAGONAL_MATRIX<T,n>& M)
{
    TENSOR<T,n,m,n> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            for(int k=0;k<n;k++)
                t.x(j)(i,k)+=a.x(i)(j,k)*M(k,k);
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==0,TENSOR<T,n,n,m> >::type
Contract(const SYMMETRIC_TENSOR<T,2,m,n>& a,const DIAGONAL_MATRIX<T,n>& M)
{
    TENSOR<T,n,n,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            for(int k=0;k<n;k++)
                t.x(j)(k,i)+=a.x(i)(j,k)*M(j,j);
    return t;
}

template<int r,int s,class T,int m,int n>
typename enable_if<r==1,TENSOR<T,n,n,m> >::type
Contract(const SYMMETRIC_TENSOR<T,2,m,n>& a,const DIAGONAL_MATRIX<T,n>& M)
{
    TENSOR<T,n,n,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            for(int k=0;k<n;k++)
                t.x(j)(k,i)+=a.x(i)(j,k)*M(k,k);
    return t;
}

template<int r,int s,class T,int m,int p>
typename enable_if<s==0,SYMMETRIC_TENSOR<T,r,p> >::type
Contract(const DIAGONAL_TENSOR<T,m>& a,const MATRIX<T,m,p>& M)
{
    SYMMETRIC_TENSOR<T,r,p,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<p;j++)
            t.x(j)(i,i)+=a.v(i)*M(i,j);
    return t;
}

template<int r,int s,class T,int m,int p>
typename enable_if<s==1,SYMMETRIC_TENSOR<T,r,p,m> >::type
Contract(const DIAGONAL_TENSOR<T,m>& a,const MATRIX<T,p,m>& M)
{
    SYMMETRIC_TENSOR<T,r,p,m> t;
    for(int i=0;i<m;i++)
        for(int j=0;j<p;j++)
            t.x(j)(i,i)+=a.v(i)*M(j,i);
    return t;
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<s==0,VEC_ID_TENSOR<T,r,p,n> >::type
Contract(const VEC_ID_TENSOR<T,r,m,n>& a,const MATRIX<T,m,p>& M)
{
    return VEC_ID_TENSOR<T,r,p,n>(M.Transpose_Times(a.v));
}

template<int r,int s,class T,int m,int n,int p>
typename enable_if<s==1,VEC_ID_TENSOR<T,r,p,n> >::type
Contract(const VEC_ID_TENSOR<T,r,m,n>& a,const MATRIX<T,p,m>& M)
{
    return VEC_ID_TENSOR<T,r,p,n>(M*a.v);
}

template<int r,int s,class T,int m,int n> VEC_ID_TENSOR<T,r,m,n>
Contract(const VEC_ID_TENSOR<T,r,m,n>& a,const SYMMETRIC_MATRIX<T,m>& M)
{
    return VEC_ID_TENSOR<T,r,m,n>(M*a.v);
}

template<int r,int s,class T,int m,int n,int p0,int p1,int u> auto
Contract(const VEC_ID_TENSOR<T,u,m,n>& a,const MATRIX<T,p0,p1>& M)
    -> typename enable_if<r!=u && (r<3-r-u)==(s==0) && n==(s==0?p0:p1),decltype(Tensor_Product<u>(M.Transposed(),a.v))>::type
{
    return Tensor_Product<u>(M.Transposed(),a.v);
}

template<int r,int s,class T,int m,int n,int p0,int p1,int u> auto
Contract(const VEC_ID_TENSOR<T,u,m,n>& a,const MATRIX<T,p0,p1>& M)
    -> typename enable_if<r!=u && (r<3-r-u)!=(s==0) && n==(s==0?p0:p1),decltype(Tensor_Product<u>(M,a.v))>::type
{
    return Tensor_Product<u>(M,a.v);
}

template<int r,int s,class T,int m,int n,int u> auto
Contract(const VEC_ID_TENSOR<T,u,m,n>& a,const SYMMETRIC_MATRIX<T,n>& M)
    -> typename enable_if<r!=u,decltype(Tensor_Product<u>(M,a.v))>::type
{
    return Tensor_Product<u>(M,a.v);
}

template<int r,int s,class T,int u,int m,int p> auto
Contract(const VEC_ID_SYM_TENSOR<T,u,m>& a,const MATRIX<T,m,p>& M)
    -> decltype(Contract<r,s>(VEC_ID_TENSOR<T,u==0?1:0,m,m>(a.v),M)+Contract<r,s>(VEC_ID_TENSOR<T,u==2?1:2,m,m>(a.v),M))
{
    return Contract<r,s>(VEC_ID_TENSOR<T,u==0?1:0,m,m>(a.v),M)+Contract<r,s>(VEC_ID_TENSOR<T,u==2?1:2,m,m>(a.v),M);
}

template<int r,int s,class T,int u,int m> auto
Contract(const VEC_ID_SYM_TENSOR<T,u,m>& a,const SYMMETRIC_MATRIX<T,m>& M)
    -> decltype(Contract<r,s>(VEC_ID_TENSOR<T,u==0?1:0,m,m>(a.v),M)+Contract<r,s>(VEC_ID_TENSOR<T,u==2?1:2,m,m>(a.v),M))
{
    return Contract<r,s>(VEC_ID_TENSOR<T,u==0?1:0,m,m>(a.v),M)+Contract<r,s>(VEC_ID_TENSOR<T,u==2?1:2,m,m>(a.v),M);
}

template<int r,int s,class T,int q>
typename enable_if<r==0 && s==0,TENSOR<T,q,3,3> >::type
Contract(const PERMUTATION_TENSOR<T>& a,const MATRIX<T,3,q>& M)
{
    TENSOR<T,q,3,3> t;
    for(int i=0;i<q;i++) t.x(i)=MATRIX<T,3>::Cross_Product_Matrix(-a.x*M.Column(i));
    return t;
}

template<int r,int s,class T,int q>
typename enable_if<r==1 && s==0,TENSOR<T,3,q,3> >::type
Contract(const PERMUTATION_TENSOR<T>& a,const MATRIX<T,3,q>& M)
{
    TENSOR<T,3,q,3> t;
    for(int i=0;i<q;i++){
        MATRIX<T,3> c=MATRIX<T,3>::Cross_Product_Matrix(a.x*M.Column(i));
        for(int j=0;j<3;j++) t.x(j).Set_Row(i,c.Row(j));}
    return t;
}

template<int r,int s,class T,int q>
typename enable_if<r==2 && s==0,TENSOR<T,3,3,q> >::type
Contract(const PERMUTATION_TENSOR<T>& a,const MATRIX<T,3,q>& M)
{
    TENSOR<T,3,3,q> t;
    for(int i=0;i<q;i++){
        MATRIX<T,3> c=MATRIX<T,3>::Cross_Product_Matrix(a.x*M.Column(i));
        for(int j=0;j<3;j++) t.x(j).Set_Column(i,c.Column(j));}
    return t;
}

template<int r,int s,class T,int q>
typename enable_if<r==0 && s==1,TENSOR<T,q,3,3> >::type
Contract(const PERMUTATION_TENSOR<T>& a,const MATRIX<T,q,3>& M)
{
    TENSOR<T,q,3,3> t;
    for(int i=0;i<q;i++) t.x(i)=MATRIX<T,3>::Cross_Product_Matrix(-a.x*M.Row(i));
    return t;
}

template<int r,int s,class T,int q>
typename enable_if<r==1 && s==1,TENSOR<T,3,q,3> >::type
Contract(const PERMUTATION_TENSOR<T>& a,const MATRIX<T,q,3>& M)
{
    TENSOR<T,3,q,3> t;
    for(int i=0;i<q;i++){
        MATRIX<T,3> c=MATRIX<T,3>::Cross_Product_Matrix(a.x*M.Row(i));
        for(int j=0;j<3;j++) t.x(j).Set_Row(i,c.Row(j));}
    return t;
}

template<int r,int s,class T,int q>
typename enable_if<r==2 && s==1,TENSOR<T,3,3,q> >::type
Contract(const PERMUTATION_TENSOR<T>& a,const MATRIX<T,q,3>& M)
{
    TENSOR<T,3,3,q> t;
    for(int i=0;i<q;i++){
        MATRIX<T,3> c=MATRIX<T,3>::Cross_Product_Matrix(a.x*M.Row(i));
        for(int j=0;j<3;j++) t.x(j).Set_Column(i,c.Column(j));}
    return t;
}

template<int r,int s,class T>
typename enable_if<r==0,TENSOR<T,3,3,3> >::type
Contract(const PERMUTATION_TENSOR<T>& a,const SYMMETRIC_MATRIX<T,3>& M)
{
    TENSOR<T,3,3,3> t;
    for(int i=0;i<3;i++) t.x(i)=MATRIX<T,3>::Cross_Product_Matrix(-a.x*M.Column(i));
    return t;
}

template<int r,int s,class T>
typename enable_if<r==1,TENSOR<T,3,3,3> >::type
Contract(const PERMUTATION_TENSOR<T>& a,const SYMMETRIC_MATRIX<T,3>& M)
{
    TENSOR<T,3,3,3> t;
    for(int i=0;i<3;i++){
        MATRIX<T,3> c=MATRIX<T,3>::Cross_Product_Matrix(a.x*M.Column(i));
        for(int j=0;j<3;j++) t.x(j).Set_Row(i,c.Row(j));}
    return t;
}

template<int r,int s,class T>
typename enable_if<r==2,TENSOR<T,3,3,3> >::type
Contract(const PERMUTATION_TENSOR<T>& a,const SYMMETRIC_MATRIX<T,3>& M)
{
    TENSOR<T,3,3,3> t;
    for(int i=0;i<3;i++){
        MATRIX<T,3> c=MATRIX<T,3>::Cross_Product_Matrix(a.x*M.Column(i));
        for(int j=0;j<3;j++) t.x(j).Set_Column(i,c.Column(j));}
    return t;
}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const ZERO_MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const ZERO_MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const ZERO_MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const ZERO_MATRIX<T,m>& cm1,const MATRIX<T,m>){return ZERO_TENSOR<T,m>();}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const IDENTITY_MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const IDENTITY_MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const IDENTITY_MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const IDENTITY_MATRIX<T,m>& cm1,const MATRIX<T,m>& cm2)
{return Contract_2(t,cm2).Twice_Symmetric_Part_12();}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const SCALE_MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const SCALE_MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const SCALE_MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const SCALE_MATRIX<T,m>& cm1,const MATRIX<T,m>& cm2)
{return Symmetric_Double_Contract_12(t*cm1.x,IDENTITY_MATRIX<T,m>(),cm2);}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2)
{return Symmetric_Double_Contract_12(-t,cm2,cm1);}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2)
{return Symmetric_Double_Contract_12(-t,cm2,cm1);}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const PERMUTATION_TENSOR<T>& t,const MATRIX<T,m>& cm1,const MATRIX<T,m>& cm2)
{
    SYMMETRIC_TENSOR<T,0,m,m> s;
    VECTOR<T,m> c10=cm1.Row(0),c11=cm1.Row(1),c12=cm1.Row(2),c20=t.x*cm2.Row(0),c21=t.x*cm2.Row(1),c22=t.x*cm2.Row(2);
    s.x(0)=SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c11,c22)-SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c12,c21);
    s.x(1)=SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c12,c20)-SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c10,c22);
    s.x(2)=SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c10,c21)-SYMMETRIC_MATRIX<T,m>::Symmetric_Outer_Product(c11,c20);
    return s;
}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const ZERO_MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const ZERO_MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const ZERO_MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const ZERO_MATRIX<T,m>& cm1,const MATRIX<T,m>){return ZERO_TENSOR<T,m>();}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const IDENTITY_MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> DIAGONAL_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const IDENTITY_MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2){return t*2;}
template<class T,int m> DIAGONAL_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const IDENTITY_MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2){return t*(cm2.x*2);}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const IDENTITY_MATRIX<T,m>& cm1,const MATRIX<T,m>& cm2)
{return (DIAGONAL_MATRIX<T,m>(t.v)*cm2).Twice_Symmetric_Part();}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const SCALE_MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> DIAGONAL_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const SCALE_MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2){return t*(cm1.x*2);}
template<class T,int m> DIAGONAL_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const SCALE_MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2){return t*(cm1.x*cm2.x*2);}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const SCALE_MATRIX<T,m>& cm1,const MATRIX<T,m>& cm2)
{return (DIAGONAL_MATRIX<T,m>(cm1.x*t.v)*cm2).Twice_Symmetric_Part();}

template<class T,int m> ZERO_TENSOR<T,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const MATRIX<T,m>& cm1,const ZERO_MATRIX<T,m>& cm2){return ZERO_TENSOR<T,m>();}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const MATRIX<T,m>& cm1,const IDENTITY_MATRIX<T,m>& cm2)
{return Symmetric_Double_Contract_12(t,cm2,cm1);}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const MATRIX<T,m>& cm1,const SCALE_MATRIX<T,m>& cm2)
{return Symmetric_Double_Contract_12(t,cm2,cm1);}
template<class T,int m> SYMMETRIC_TENSOR<T,0,m,m> Symmetric_Double_Contract_12(const DIAGONAL_TENSOR<T,m>& t,const MATRIX<T,m>& cm1,const MATRIX<T,m>& cm2)
{
    SYMMETRIC_TENSOR<T,0,m,m> s;
    for(int i=0;i<m;i++) s.x(i)=t.v(i)*MATRIX<T,m>::Outer_Product(cm1.Row(i),cm2.Row(i)).Twice_Symmetric_Part();;
    return s;
}

template<class T_TEN> typename enable_if<IS_TENSOR<T_TEN>::value,T_TEN>::type Choose(const T_TEN& a,const T_TEN& b);
template<class T_TEN0,class T_TEN1>
typename enable_if<IS_TENSOR<T_TEN0>::value &&
    IS_TENSOR<T_TEN1>::value &&
    !(IS_SYM_TENSOR<0,T_TEN0>::value && IS_SYM_TENSOR<0,T_TEN1>::value) &&
    !(IS_SYM_TENSOR<1,T_TEN0>::value && IS_SYM_TENSOR<1,T_TEN1>::value) &&
    !(IS_SYM_TENSOR<2,T_TEN0>::value && IS_SYM_TENSOR<2,T_TEN1>::value),
    TENSOR<typename T_TEN0::SCALAR,T_TEN0::m,T_TEN0::n,T_TEN0::p> >::TYPE
Choose(const T_TEN0& a,const T_TEN1& b);
template<class T_TEN0,class T_TEN1> typename enable_if<IS_SYM_TENSOR<0,T_TEN0>::value && IS_SYM_TENSOR<0,T_TEN1>::value,SYMMETRIC_TENSOR<typename T_TEN0::SCALAR,0,T_TEN0::um,T_TEN0::un> >::type Choose(const T_TEN0& a,const T_TEN1& b);
template<class T_TEN0,class T_TEN1> typename enable_if<IS_SYM_TENSOR<1,T_TEN0>::value && IS_SYM_TENSOR<1,T_TEN1>::value,SYMMETRIC_TENSOR<typename T_TEN0::SCALAR,1,T_TEN0::um,T_TEN0::un> >::type Choose(const T_TEN0& a,const T_TEN1& b);
template<class T_TEN0,class T_TEN1> typename enable_if<IS_SYM_TENSOR<2,T_TEN0>::value && IS_SYM_TENSOR<2,T_TEN1>::value,SYMMETRIC_TENSOR<typename T_TEN0::SCALAR,2,T_TEN0::um,T_TEN0::un> >::type Choose(const T_TEN0& a,const T_TEN1& b);
template<class T,int m,class T_TEN> T_TEN Choose(const T_TEN& a,const ZERO_TENSOR<T,T_TEN::m,T_TEN::n,T_TEN::p>& b);
template<class T,int m,class T_TEN> T_TEN Choose(const ZERO_TENSOR<T,T_TEN::m,T_TEN::n,T_TEN::p>& a,const T_TEN& b);
template<class T,int m,int n,int p> ZERO_TENSOR<T,m,n,p> Choose(const ZERO_TENSOR<T,m,n,p>& a,const ZERO_TENSOR<T,m,n,p>& b);

template<class T> typename enable_if<IS_TENSOR<T>::value||IS_MATRIX<T>::value||IS_VECTOR<T>::value||is_scalar<T>::value,T>::type Choose_Zero(const T& a);
template<class T,int d> SCALE_MATRIX<T,d> Choose_Zero(const IDENTITY_MATRIX<T,d>& a);
template<class T,int m> typename conditional<m==0,FIXED_NUMBER<T,m>,T>::type Choose_Zero(const FIXED_NUMBER<T,m>& a);

//#####################################################################
// operator=
//##################################################################### 

template<class T> typename enable_if<IS_TENSOR<T>::value>::type Fill_From(T& a,const T& b){a=b;}

template<class T,int m,int n,int p> void
Fill_From(TENSOR<T,m,n,p>& a,const ZERO_TENSOR<T,m,n,p>& b)
{
    a=TENSOR<T,m,n,p>();
}

template<class T,int m,int n> void
Fill_From(TENSOR<T,m,n,n>& a,const SYMMETRIC_TENSOR<T,0,m,n>& b)
{
    for(int i=0;i<m;i++) a.x(i)=b.x(i);
}

template<class T,int m,int n> void
Fill_From(TENSOR<T,n,m,n>& a,const SYMMETRIC_TENSOR<T,1,m,n>& b)
{
    a=TENSOR<T,n,m,n>();
    a+=b;
}

template<class T,int m,int n> void
Fill_From(TENSOR<T,n,n,m>& a,const SYMMETRIC_TENSOR<T,2,m,n>& b)
{
    a=TENSOR<T,n,n,m>();
    a+=b;
}

template<class T,int m,int n> void
Fill_From(TENSOR<T,m,n,n>& a,const VEC_ID_TENSOR<T,0,m,n>& b)
{
    a=TENSOR<T,m,n,n>();
    a+=b;
}

template<class T,int m,int n> void
Fill_From(TENSOR<T,n,m,n>& a,const VEC_ID_TENSOR<T,1,m,n>& b)
{
    a=TENSOR<T,n,m,n>();
    a+=b;
}

template<class T,int m,int n> void
Fill_From(TENSOR<T,n,n,m>& a,const VEC_ID_TENSOR<T,2,m,n>& b)
{
    a=TENSOR<T,n,n,m>();
    a+=b;
}

template<class T,int u,int m> void
Fill_From(TENSOR<T,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    a=TENSOR<T,m>();
    a+=b;
}

template<class T> void
Fill_From(TENSOR<T,3,3,3>& a,const PERMUTATION_TENSOR<T>& b)
{
    a=TENSOR<T,3,3,3>();
    a+=b;
}

template<class T,int u,int m,int n> void
Fill_From(SYMMETRIC_TENSOR<T,u,m,n>& a,const ZERO_TENSOR<T,u==0?m:n,u==1?m:n,u==2?m:n>& b)
{
    a=SYMMETRIC_TENSOR<T,u,m,n>();
}

template<class T,int u,int m,int n> void
Fill_From(SYMMETRIC_TENSOR<T,u,m,n>& a,const VEC_ID_TENSOR<T,u,m,n>& b)
{
    a=SYMMETRIC_TENSOR<T,u,m,n>();
    a+=b;
}

template<class T,int u,int m> void
Fill_From(SYMMETRIC_TENSOR<T,u,m>& a,const VEC_ID_SYM_TENSOR<T,u,m>& b)
{
    TENSOR<T,m> t;
    t+=VEC_ID_TENSOR<T,u==0?1:0,m,m>(b.v);
    for(int i=0;i<m;i++) a.x(i)=t.x(i).Twice_Symmetric_Part();
}

template<class T,int u,int m> void
Fill_From(SYMMETRIC_TENSOR<T,u,m,m>& a,const DIAGONAL_TENSOR<T,m>& b)
{
    for(int i=0;i<m;i++) a.x(i)(i,i)=b.v(i);
}
}
#endif

