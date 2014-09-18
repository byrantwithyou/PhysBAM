//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRIMITIVE_MATRICES
//##################################################################### 
#ifndef __PRIMITIVE_MATRICES__
#define __PRIMITIVE_MATRICES__

#include <Tools/Matrices/IDENTITY_MATRIX.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Vectors/VECTOR.h>
#include <Tools/Vectors/ZERO_VECTOR.h>
#include <cmath>
namespace PhysBAM{
template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(ZERO_MATRIX<T,d> z){return SYMMETRIC_MATRIX<T,d>();}
template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(SCALE_MATRIX<T,d> z){return SYMMETRIC_MATRIX<T,d>()+z.x;}
template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(IDENTITY_MATRIX<T,d> z){return SYMMETRIC_MATRIX<T,d>::Identity_Matrix();}
template<class T,int d> inline SYMMETRIC_MATRIX<T,d> Cast_Helper(SYMMETRIC_MATRIX<T,d> z){return z;}
template<class T,int d> inline MATRIX<T,d> Cast_Helper(const MATRIX<T,d>& z){return z;}

template<class T> typename ENABLE_IF<IS_MATRIX<T>::value>::TYPE Fill_From(T& a,const T& b){a=b;}
template<class T,int d> void Fill_From(MATRIX<T,d>& m,const SYMMETRIC_MATRIX<T,d>& z){m=z;}
template<class T,int d> void Fill_From(MATRIX<T,d>& m,const ZERO_MATRIX<T,d>& z){m=MATRIX<T,d>();}
template<class T,int d> void Fill_From(MATRIX<T,d>& m,const IDENTITY_MATRIX<T,d>& i){m=MATRIX<T,d>::Identity_Matrix();}
template<class T,int d> void Fill_From(MATRIX<T,d>& m,const SCALE_MATRIX<T,d>& s){m=MATRIX<T,d>()+s.x;}
template<class T,int d> void Fill_From(SYMMETRIC_MATRIX<T,d>& m,const ZERO_MATRIX<T,d>& z){m=SYMMETRIC_MATRIX<T,d>();}
template<class T,int d> void Fill_From(SYMMETRIC_MATRIX<T,d>& m,const IDENTITY_MATRIX<T,d>& i){m=SYMMETRIC_MATRIX<T,d>::Identity_Matrix();}
template<class T,int d> void Fill_From(SYMMETRIC_MATRIX<T,d>& m,const SCALE_MATRIX<T,d>& s){m=SYMMETRIC_MATRIX<T,d>()+s.x;}
template<class T,int d> void Fill_From(SCALE_MATRIX<T,d>& s,const ZERO_MATRIX<T,d>& z) {s.x=0;}
template<class T,int d> void Fill_From(SCALE_MATRIX<T,d>& s,const IDENTITY_MATRIX<T,d>& i) {s.x=1;}

template<class T_MAT> typename ENABLE_IF<IS_MATRIX<T_MAT>::value,T_MAT>::TYPE Choose(const T_MAT& a,const T_MAT& b);
template<class T_MAT,class T_MAT1> typename ENABLE_IF<IS_MATRIX<T_MAT>::value&&IS_MATRIX<T_MAT1>::value&&(!IS_SYMMETRIC_MATRIX<T_MAT>::value || !IS_SYMMETRIC_MATRIX<T_MAT1>::value),MATRIX<typename T_MAT::SCALAR,T_MAT::m> >::TYPE Choose(const T_MAT& a,const T_MAT1& b);
template<class T_MAT,class T_MAT1> typename ENABLE_IF<IS_SYMMETRIC_MATRIX<T_MAT>::value && IS_SYMMETRIC_MATRIX<T_MAT1>::value,SYMMETRIC_MATRIX<typename T_MAT::SCALAR,T_MAT::m> >::TYPE Choose(const T_MAT& a,const T_MAT1& b);
template<class T,int d,class T_MAT> T_MAT Choose(const T_MAT& a,const ZERO_MATRIX<T,d>& b);
template<class T,int d,class T_MAT> T_MAT Choose(const ZERO_MATRIX<T,d>& a,const T_MAT& b);
template<class T,int d> ZERO_MATRIX<T,d> Choose(const ZERO_MATRIX<T,d>& a,const ZERO_MATRIX<T,d>& b);
template<class T,int d> SCALE_MATRIX<T,d> Choose(const SCALE_MATRIX<T,d>& a,const IDENTITY_MATRIX<T,d>& b);
template<class T,int d> SCALE_MATRIX<T,d> Choose(const IDENTITY_MATRIX<T,d>& a,const SCALE_MATRIX<T,d>& b);
}
#endif
