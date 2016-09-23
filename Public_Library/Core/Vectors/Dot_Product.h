//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Dot_Product 
//#####################################################################
#ifndef __Dot_Product__
#define __Dot_Product__

#include <Core/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{

template<class TV>
inline typename TV::SCALAR Dot_Product(const TV& v1,const TV& v2)
{return v1.Dot(v2);}

template<class T>
inline typename enable_if<is_scalar<T>::value,T>::type Dot_Product(const T a1,const T a2)
{return a1*a2;}

template<class TV>
inline double Dot_Product_Double_Precision(const TV& v1,const TV& v2)
{return v1.Dot_Double_Precision(v2);}

template<class T,int d> class VECTOR;
template<class T,int d>
inline double Dot_Product_Double_Precision(const VECTOR<T,d>& v1,const VECTOR<T,d>& v2)
{return v1.Dot(v2);}

template<class T>
inline typename enable_if<is_scalar<T>::value,double>::type Dot_Product_Double_Precision(const T a1,const T a2)
{return a1*a2;}

template<class T>
inline typename enable_if<is_scalar<T>::value,T>::type Magnitude_Squared(const T a)
{return a*a;}

template<class TV>
inline typename TV::SCALAR Magnitude_Squared(const TV& v)
{return v.Magnitude_Squared();}
}
#endif
