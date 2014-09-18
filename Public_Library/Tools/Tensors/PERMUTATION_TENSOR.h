//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PERMUTATION_TENSOR
//##################################################################### 
#ifndef __PERMUTATION_TENSOR__
#define __PERMUTATION_TENSOR__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// T_ijk=x*e_ijk
template<class T>
struct PERMUTATION_TENSOR
{
    typedef T SCALAR;
    enum {m=3,n=3,p=3};
    T x;

    explicit PERMUTATION_TENSOR(T x=T()): x(x) {}

    PERMUTATION_TENSOR operator-() const
    {return PERMUTATION_TENSOR(-x);}

    PERMUTATION_TENSOR operator*(T a) const
    {return PERMUTATION_TENSOR(x*a);}

    PERMUTATION_TENSOR operator/(T a) const
    {return PERMUTATION_TENSOR(x/a);}
};
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
