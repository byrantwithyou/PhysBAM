//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONAL_TENSOR
//##################################################################### 
#ifndef __DIAGONAL_TENSOR__
#define __DIAGONAL_TENSOR__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// T_iii = u_i, T_ijk = 0 otherwise
template<class T,int mm>
class DIAGONAL_TENSOR
{
public:
    static const bool is_tensor=true;
    typedef T SCALAR;
    enum {m=mm,n=mm,p=mm,um=mm,un=mm};

    VECTOR<T,m> v;
    explicit DIAGONAL_TENSOR(VECTOR<T,m> v=(VECTOR<T,m>())): v(v) {}

    DIAGONAL_TENSOR operator-() const
    {return DIAGONAL_TENSOR(-v);}

    DIAGONAL_TENSOR operator*(T a) const
    {return DIAGONAL_TENSOR(v*a);}

    DIAGONAL_TENSOR operator/(T a) const
    {return DIAGONAL_TENSOR(v/a);}
};
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
