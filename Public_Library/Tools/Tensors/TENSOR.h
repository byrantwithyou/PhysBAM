//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TENSOR
//##################################################################### 
#ifndef __TENSOR__
#define __TENSOR__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/SYMMETRIC_TENSOR.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// general T_ijk
template<class T,int mm,int nn,int pp>
struct TENSOR
{
    typedef T SCALAR;
    enum {m=mm,n=nn,p=pp};
    VECTOR<MATRIX<T,n,p>,m> x;

    TENSOR(){}

    TENSOR(const SYMMETRIC_TENSOR<T,m,n>& t)
    {for(int i=0;i<m;i++) x(i)=t.x(i);}

    TENSOR operator-() const
    {TENSOR t;for(int i=0;i<m;i++) t.x(i)=-x(i);return t;}

    TENSOR operator*(T a) const
    {TENSOR t;for(int i=0;i<m;i++) t.x(i)=x(i)*a;return t;}

    TENSOR operator/(T a) const
    {TENSOR t;for(int i=0;i<m;i++) t.x(i)=x(i)/a;return t;}

    SYMMETRIC_TENSOR<T,m,n> Twice_Symmetric_Part_12() const
    {SYMMETRIC_TENSOR<T,m,n> t;for(int i=0;i<m;i++) t.x(i)=x(i).Twice_Symmetric_Part();return t;}

    TENSOR<T,m,p,n> Transposed() const // swaps last two indices
    {TENSOR<T,m,p,n> t;for(int i=0;i<m;i++) t.x(i)=x(i).Transposed();return t;}
};
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif