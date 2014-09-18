//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TENSORS_FORWARD
//##################################################################### 
#ifndef __TENSORS_FORWARD__
#define __TENSORS_FORWARD__

#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Tensors/TENSORS_FORWARD.h>
#include <Tools/Vectors/VECTOR.h>
#include <cmath>
namespace PhysBAM{

// general T_ijk
template<class T,int mm,int nn=mm,int pp=nn> struct TENSOR;
template<class T,int mm,int nn=mm> struct SYMMETRIC_TENSOR;
template<class T,int mm,int nn=mm,int pp=nn> struct ZERO_TENSOR;
template<class T,int mm,int nn=mm> struct VEC_ID_TENSOR_0;
template<class T,int mm,int nn=mm> struct VEC_ID_TENSOR_1;
template<class T,int mm,int nn=mm> struct VEC_ID_TENSOR_2;
template<class T,int mm> struct VEC_ID_TENSOR_12;
template<class T> struct PERMUTATION_TENSOR;

template<class T> struct IS_TENSOR{static const int value=0;};
template<class T,int m,int n,int p> struct IS_TENSOR<ZERO_TENSOR<T,m,n,p> > {static const int value=1;};
template<class T,int m,int n,int p> struct IS_TENSOR<TENSOR<T,m,n,p> > {static const int value=1;};
template<class T,int m,int n> struct IS_TENSOR<SYMMETRIC_TENSOR<T,m,n> > {static const int value=1;};
template<class T,int m,int n> struct IS_TENSOR<VEC_ID_TENSOR_0<T,m,n> > {static const int value=1;};
template<class T,int m,int n> struct IS_TENSOR<VEC_ID_TENSOR_1<T,m,n> > {static const int value=1;};
template<class T,int m,int n> struct IS_TENSOR<VEC_ID_TENSOR_2<T,m,n> > {static const int value=1;};
template<class T,int m> struct IS_TENSOR<VEC_ID_TENSOR_12<T,m> > {static const int value=1;};
template<class T> struct IS_TENSOR<PERMUTATION_TENSOR<T> > {static const int value=1;};

template<class T> struct IS_SYM_TENSOR{static const int value=0;};
template<class T,int m,int n> struct IS_SYM_TENSOR<ZERO_TENSOR<T,m,n,n> > {static const int value=1;};
template<class T,int m,int n> struct IS_SYM_TENSOR<SYMMETRIC_TENSOR<T,m,n> > {static const int value=1;};
template<class T,int m,int n> struct IS_SYM_TENSOR<VEC_ID_TENSOR_0<T,m,n> > {static const int value=1;};
template<class T,int m> struct IS_SYM_TENSOR<VEC_ID_TENSOR_12<T,m> > {static const int value=1;};
}
#endif
