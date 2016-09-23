//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ZERO_TENSOR
//##################################################################### 
#ifndef __ZERO_TENSOR__
#define __ZERO_TENSOR__

#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <cmath>
namespace PhysBAM{

// T_ijk=0
template<class T,int mm,int nn,int pp>
class ZERO_TENSOR
{
public:
    static const bool is_tensor=true;
    typedef T SCALAR;
    enum {m=mm,n=nn,p=pp};
    ZERO_TENSOR operator-() const
    {return *this;}

    ZERO_TENSOR operator*(T a) const
    {return *this;}

    ZERO_TENSOR operator/(T a) const
    {return *this;}
};
template<class T,int m,int n,int p> inline std::ostream&
operator<<(std::ostream& o,const ZERO_TENSOR<T,m,n,p>& A)
{o<<"(0)";return o;}
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
