//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_TENSOR
//##################################################################### 
#ifndef __SYMMETRIC_TENSOR__
#define __SYMMETRIC_TENSOR__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// symmetric T_ijk=T_ikj
template<class T,int mm,int nn>
struct SYMMETRIC_TENSOR
{
    typedef T SCALAR;
    enum {m=mm,n=nn,p=nn};
    VECTOR<SYMMETRIC_MATRIX<T,n>,m> x;

    SYMMETRIC_TENSOR operator-() const
    {SYMMETRIC_TENSOR t;for(int i=0;i<m;i++) t.x(i)=-x(i);return t;}

    SYMMETRIC_TENSOR operator*(T a) const
    {SYMMETRIC_TENSOR t;for(int i=0;i<m;i++) t.x(i)=x(i)*a;return t;}

    SYMMETRIC_TENSOR operator/(T a) const
    {SYMMETRIC_TENSOR t;for(int i=0;i<m;i++) t.x(i)=x(i)/a;return t;}
};
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
