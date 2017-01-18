//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TENSOR
//##################################################################### 
#ifndef __TENSOR__
#define __TENSOR__

#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <cmath>
namespace PhysBAM{

// general T_ijk
template<class T,int mm,int nn,int pp>
class TENSOR
{
public:
    static const bool is_tensor=true;
    typedef T SCALAR;
    enum {m=mm,n=nn,p=pp};
    VECTOR<MATRIX<T,n,p>,m> x;

    TENSOR(){}

    const T& operator()(int i,int j,int k) const
    {return x(i)(j,k);}

    T& operator()(int i,int j,int k)
    {return x(i)(j,k);}

    TENSOR operator-() const
    {TENSOR t;
    for(int i=0;i<m;i++) t.x(i)=-x(i);
    return t;}

    TENSOR operator*(T a) const
    {TENSOR t;
    for(int i=0;i<m;i++) t.x(i)=x(i)*a;
    return t;}

    TENSOR operator/(T a) const
    {TENSOR t;
    for(int i=0;i<m;i++) t.x(i)=x(i)/a;
    return t;}

    SYMMETRIC_TENSOR<T,m,n> Twice_Symmetric_Part_12() const
    {SYMMETRIC_TENSOR<T,m,n> t;
    for(int i=0;i<m;i++) t.x(i)=x(i).Twice_Symmetric_Part();
    return t;}

    TENSOR<T,m,p,n> Transposed() const // swaps last two indices
    {TENSOR<T,m,p,n> t;
    for(int i=0;i<m;i++) t.x(i)=x(i).Transposed();
    return t;}
};

template<class T,int m,int n,int p>
TENSOR<T,m,n,p> operator*(T a,const TENSOR<T,m,n,p>& s)
{return s*a;}

template<class T,int m,int n,int p> inline std::ostream&
operator<<(std::ostream& o,const TENSOR<T,m,n,p>& A)
{o<<"("<<A.x<<")";return o;}
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
