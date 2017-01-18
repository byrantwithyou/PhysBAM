//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_TENSOR
//##################################################################### 
#ifndef __SYMMETRIC_TENSOR__
#define __SYMMETRIC_TENSOR__

#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <cmath>
namespace PhysBAM{

// uu = unique index (0..2)
// uu=0: T_ijk=T_ikj
// uu=1: T_ijk=T_kji
// uu=2: T_ijk=T_jik
template<class T,int uu,int mm,int nn>
class SYMMETRIC_TENSOR
{
public:
    static const bool is_tensor=true;
    STATIC_ASSERT(uu>=0 && uu<=2);
    typedef T SCALAR;
    enum {u=uu,m=uu==0?mm:nn,n=uu==1?mm:nn,p=uu==2?mm:nn,um=mm,un=nn};
    VECTOR<SYMMETRIC_MATRIX<T,nn>,mm> x;

    SYMMETRIC_TENSOR operator-() const
    {SYMMETRIC_TENSOR t;
    for(int i=0;i<mm;i++) t.x(i)=-x(i);
    return t;}

    SYMMETRIC_TENSOR operator*(T a) const
    {SYMMETRIC_TENSOR t;
    for(int i=0;i<mm;i++) t.x(i)=x(i)*a;
    return t;}

    SYMMETRIC_TENSOR operator/(T a) const
    {SYMMETRIC_TENSOR t;
    for(int i=0;i<mm;i++) t.x(i)=x(i)/a;
    return t;}
};

template<class T,int u,int m,int n>
SYMMETRIC_TENSOR<T,u,m,n> operator*(T a,const SYMMETRIC_TENSOR<T,u,m,n>& s)
{return s*a;}

template<class T,int m,int n,int p> inline std::ostream&
operator<<(std::ostream& o,const SYMMETRIC_TENSOR<T,m,n,p>& A)
{o<<"("<<A.x<<")";return o;}
}
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#endif
