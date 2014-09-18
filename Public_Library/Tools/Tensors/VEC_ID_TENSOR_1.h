//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VEC_ID_TENSOR_1
//##################################################################### 
#ifndef __VEC_ID_TENSOR_1__
#define __VEC_ID_TENSOR_1__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// T_ijk=v_j*delta_ik
template<class T,int mm,int nn>
struct VEC_ID_TENSOR_1
{
    typedef T SCALAR;
    enum {m=mm,n=nn,p=mm};

    VECTOR<T,n> v;
    explicit VEC_ID_TENSOR_1(VECTOR<T,n> v=VECTOR<T,n>()): v(v) {}

    VEC_ID_TENSOR_1 operator-() const
    {return VEC_ID_TENSOR_1(-v);}

    VEC_ID_TENSOR_1 operator*(T a) const
    {return VEC_ID_TENSOR_1(v*a);}

    VEC_ID_TENSOR_1 operator/(T a) const
    {return VEC_ID_TENSOR_1(v/a);}
};
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
