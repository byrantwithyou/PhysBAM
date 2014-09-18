//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VEC_ID_TENSOR_12
//##################################################################### 
#ifndef __VEC_ID_TENSOR_12__
#define __VEC_ID_TENSOR_12__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// T_ijk=v_j*delta_ik+v_k*delta_ij
template<class T,int mm>
struct VEC_ID_TENSOR_12
{
    typedef T SCALAR;
    enum {m=mm,n=mm,p=mm};

    VECTOR<T,m> v;
    explicit VEC_ID_TENSOR_12(VECTOR<T,m> v=VECTOR<T,m>()): v(v) {}

    VEC_ID_TENSOR_12 operator-() const
    {return VEC_ID_TENSOR_12(-v);}

    VEC_ID_TENSOR_12 operator*(T a) const
    {return VEC_ID_TENSOR_12(v*a);}

    VEC_ID_TENSOR_12 operator/(T a) const
    {return VEC_ID_TENSOR_12(v/a);}
};
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
