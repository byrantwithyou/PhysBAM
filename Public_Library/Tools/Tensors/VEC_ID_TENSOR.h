//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VEC_ID_TENSOR
//##################################################################### 
#ifndef __VEC_ID_TENSOR__
#define __VEC_ID_TENSOR__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// uu = unique index (0..2)
// uu = 0: T_ijk=v_i delta_jk
// uu = 1: T_ijk=v_j delta_ik
// uu = 2: T_ijk=v_k delta_ij
template<class T,int uu,int mm,int nn>
class VEC_ID_TENSOR
{
public:
    static const bool is_tensor=true;
    STATIC_ASSERT(uu>=0 && uu<=2);
    typedef T SCALAR;
    enum {u=uu,m=uu==0?mm:nn,n=uu==1?mm:nn,p=uu==2?mm:nn,um=mm,un=nn};

    VECTOR<T,mm> v;
    explicit VEC_ID_TENSOR(VECTOR<T,mm> v=(VECTOR<T,mm>())): v(v) {}

    VEC_ID_TENSOR operator-() const
    {return VEC_ID_TENSOR(-v);}

    VEC_ID_TENSOR operator*(T a) const
    {return VEC_ID_TENSOR(v*a);}

    VEC_ID_TENSOR operator/(T a) const
    {return VEC_ID_TENSOR(v/a);}
};
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
