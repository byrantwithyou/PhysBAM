//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ZERO_TENSOR
//##################################################################### 
#ifndef __ZERO_TENSOR__
#define __ZERO_TENSOR__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// T_ijk=0
template<class T,int mm,int nn,int pp>
struct ZERO_TENSOR
{
    typedef T SCALAR;
    enum {m=mm,n=nn,p=pp};
    ZERO_TENSOR operator-() const
    {return *this;}

    ZERO_TENSOR operator*(T a) const
    {return *this;}

    ZERO_TENSOR operator/(T a) const
    {return *this;}
};
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
