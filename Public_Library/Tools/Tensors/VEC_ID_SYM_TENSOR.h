//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VEC_ID_SYM_TENSOR
//##################################################################### 
#ifndef __VEC_ID_SYM_TENSOR__
#define __VEC_ID_SYM_TENSOR__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// uu = 0: T_ijk=v_j delta_ik+v_k delta_ij
// uu = 1: T_ijk=v_i delta_jk+v_k delta_ij
// uu = 2: T_ijk=v_i delta_jk+v_j delta_ik
template<class T,int uu,int mm>
class VEC_ID_SYM_TENSOR
{
public:
    static const bool is_tensor=true;
    typedef T SCALAR;
    enum {u=uu,m=mm,n=mm,p=mm,um=mm,un=mm};

    VECTOR<T,m> v;
    VEC_ID_SYM_TENSOR() {}

    explicit VEC_ID_SYM_TENSOR(const VECTOR<T,m>& v): v(v) {}

    VEC_ID_SYM_TENSOR operator-() const
    {return VEC_ID_SYM_TENSOR(-v);}

    VEC_ID_SYM_TENSOR operator*(T a) const
    {return VEC_ID_SYM_TENSOR(v*a);}

    VEC_ID_SYM_TENSOR operator/(T a) const
    {return VEC_ID_SYM_TENSOR(v/a);}
};
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
