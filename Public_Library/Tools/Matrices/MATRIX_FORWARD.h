//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header MATRIX_FORWARD
//#####################################################################
#ifndef __MATRIX_FORWARD__
#define __MATRIX_FORWARD__

#include <Tools/Utilities/STATIC_ASSERT.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class T,class T_MATRIX> class MATRIX_BASE;
template<class T,int m_input,int n_input=m_input> class MATRIX;
template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d> class SYMMETRIC_MATRIX;
template<class T,int d> class UPPER_TRIANGULAR_MATRIX;
template<class T> class MATRIX_MXN;
template<class T> class SYMMETRIC_MATRIX_NXN;
template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class T_MATRIX> class TRANSPOSE_MATRIX;
template<class T,int mm,int nn=mm> class ZERO_MATRIX;
template<class T,int d> class IDENTITY_MATRIX;
template<class T,int d> class SCALE_MATRIX;

template<class T> struct is_scalar_BLOCK;
template<class T> struct is_scalar_VECTOR_SPACE;
template<class T,int m,int n> struct is_scalar_BLOCK<MATRIX<T,m,n> > {static const bool value=is_scalar_BLOCK<T>::value;};
template<class T,int m,int n> struct is_scalar_VECTOR_SPACE<MATRIX<T,m,n> > {static const bool value=is_scalar_VECTOR_SPACE<T>::value;};
template<class T,int m,int n,class RW> struct IS_BINARY_IO_SAFE<MATRIX<T,m,n>,RW> {static const bool value=(m>0) && (n>0) && IS_BINARY_IO_SAFE<T,RW>::value;};
template<class T,int m,int n,class SCALAR> struct REPLACE_FLOATING_POINT<MATRIX<T,m,n>,SCALAR> {typedef MATRIX<typename REPLACE_FLOATING_POINT<T,SCALAR>::TYPE,m,n> TYPE;};
template<class T,int d,class SCALAR> struct REPLACE_FLOATING_POINT<DIAGONAL_MATRIX<T,d>,SCALAR> {typedef DIAGONAL_MATRIX<typename REPLACE_FLOATING_POINT<T,SCALAR>::TYPE,d> TYPE;};
template<class T,int d,class SCALAR> struct REPLACE_FLOATING_POINT<SYMMETRIC_MATRIX<T,d>,SCALAR> {typedef SYMMETRIC_MATRIX<typename REPLACE_FLOATING_POINT<T,SCALAR>::TYPE,d> TYPE;};

template<class T_MATRIX> struct IS_MATRIX {static const bool value=false;};
template<class T> struct IS_MATRIX<MATRIX_MXN<T> > {static const bool value=true;};
template<class T> struct IS_MATRIX<SPARSE_MATRIX_FLAT_MXN<T> > {static const bool value=true;};
template<class T,int m,int n> struct IS_MATRIX<ZERO_MATRIX<T,m,n> > {static const bool value=true;};
template<class T,int m,int n> struct IS_MATRIX<MATRIX<T,m,n> > {static const bool value=true;};
template<class T,int d> struct IS_MATRIX<SYMMETRIC_MATRIX<T,d> > {static const bool value=true;};
template<class T,int d> struct IS_MATRIX<IDENTITY_MATRIX<T,d> > {static const bool value=true;};
template<class T,int d> struct IS_MATRIX<SCALE_MATRIX<T,d> > {static const bool value=true;};

template<class T_MATRIX> struct IS_SYMMETRIC_MATRIX {static const bool value=false;};
template<class T,int d> struct IS_SYMMETRIC_MATRIX<ZERO_MATRIX<T,d,d> > {static const bool value=true;};
template<class T,int d> struct IS_SYMMETRIC_MATRIX<SYMMETRIC_MATRIX<T,d> > {static const bool value=true;};
template<class T,int d> struct IS_SYMMETRIC_MATRIX<IDENTITY_MATRIX<T,d> > {static const bool value=true;};
template<class T,int d> struct IS_SYMMETRIC_MATRIX<SCALE_MATRIX<T,d> > {static const bool value=true;};

}
#endif
