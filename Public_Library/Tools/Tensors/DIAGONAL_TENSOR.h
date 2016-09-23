//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONAL_TENSOR
//##################################################################### 
#ifndef __DIAGONAL_TENSOR__
#define __DIAGONAL_TENSOR__

#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
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

template<class T,int m>
DIAGONAL_TENSOR<T,m> operator*(T a,const DIAGONAL_TENSOR<T,m>& s)
{return s*a;}

template<class T,int d> inline std::ostream&
operator<<(std::ostream& o,const DIAGONAL_TENSOR<T,d>& A)
{o<<"("<<A.v<<")";return o;}
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
