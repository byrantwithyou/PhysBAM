//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TENSORS_FORWARD
//##################################################################### 
#ifndef __TENSORS_FORWARD__
#define __TENSORS_FORWARD__

#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <cmath>
namespace PhysBAM{

// general T_ijk
template<class T,int mm,int nn=mm,int pp=nn> class TENSOR;
template<class T,int uu,int mm,int nn=mm> class SYMMETRIC_TENSOR;
template<class T,int mm,int nn=mm,int pp=nn> class ZERO_TENSOR;
template<class T,int uu,int mm,int nn=mm> class VEC_ID_TENSOR;
template<class T,int uu,int mm> class VEC_ID_SYM_TENSOR;
template<class T> class PERMUTATION_TENSOR;
template<class T,int m> class DIAGONAL_TENSOR;

template<class T> struct IS_TENSOR{static const int value=0;};
template<class T,int m,int n,int p> struct IS_TENSOR<ZERO_TENSOR<T,m,n,p> > {static const int value=1;};
template<class T,int m,int n,int p> struct IS_TENSOR<TENSOR<T,m,n,p> > {static const int value=1;};
template<class T,int u,int m,int n> struct IS_TENSOR<SYMMETRIC_TENSOR<T,u,m,n> > {static const int value=1;};
template<class T,int u,int m,int n> struct IS_TENSOR<VEC_ID_TENSOR<T,u,m,n> > {static const int value=1;};
template<class T,int u,int m> struct IS_TENSOR<VEC_ID_SYM_TENSOR<T,u,m> > {static const int value=1;};
template<class T> struct IS_TENSOR<PERMUTATION_TENSOR<T> > {static const int value=1;};
template<class T,int m> struct IS_TENSOR<DIAGONAL_TENSOR<T,m> > {static const int value=1;};

template<int s,class T> struct IS_SYM_TENSOR{static const int value=0;};
template<class T,int m,int n> struct IS_SYM_TENSOR<0,ZERO_TENSOR<T,m,n,n> > {static const int value=1;};
template<class T,int m,int n> struct IS_SYM_TENSOR<1,ZERO_TENSOR<T,n,m,n> > {static const int value=1;};
template<class T,int m,int n> struct IS_SYM_TENSOR<2,ZERO_TENSOR<T,n,n,m> > {static const int value=1;};
template<int s,class T,int m,int n> struct IS_SYM_TENSOR<s,SYMMETRIC_TENSOR<T,s,m,n> > {static const int value=1;};
template<int s,class T,int m,int n> struct IS_SYM_TENSOR<s,VEC_ID_TENSOR<T,s,m,n> > {static const int value=1;};
template<int s,class T,int m> struct IS_SYM_TENSOR<s,VEC_ID_SYM_TENSOR<T,s,m> > {static const int value=1;};
template<class T,int m,int s> struct IS_SYM_TENSOR<s,DIAGONAL_TENSOR<T,m> > {static const int value=1;};
}
#endif
