//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRAPEZOID_INTERSECTION
//##################################################################### 
#ifndef __TRAPEZOID_INTERSECTION__
#define __TRAPEZOID_INTERSECTION__    

#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
template<class T,class TV> T Trapezoid_Intersection_Area_Case_1ou(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area_Case_1oo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area_Case_1uu(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area_Case_1uo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area_Case_2ou(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area_Case_2oo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area_Case_2uu(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area_Case_2uo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area_Case_1(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area_Case_2(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
template<class T,class TV> T Trapezoid_Intersection_Area(TV a,TV b,TV c,TV d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H);
}
#endif

