//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_ORIGIN_AREAS
//##################################################################### 
#ifndef __SEGMENT_ORIGIN_AREAS__
#define __SEGMENT_ORIGIN_AREAS__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
namespace SEGMENT_ORIGIN_AREAS
{
enum POINT_CASE {inside, beyond, outside};
template<class TV> POINT_CASE Classify_Point(const TV& A,const TV& B,const TV& P);

template<class T,int n>
struct DATA
{
    T V;
    T G[2*n];
    T H[2*n][2*n];
};

template<class T,int n> void Clear(DATA<T,n>& data);

template<class TV> void Data_From_Dof(DATA<TV,1>& data,const TV& A);

// Intersect AB with OP
template<class TV> void Intersect_Segment_Point(DATA<TV,3>& data,const TV& A,const TV& B,const TV& P);

// Intersect AB with CD
template<class TV> void Intersect_Segments(DATA<TV,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);

// Area of OAB
template<class T,class TV> void Area_From_Points(DATA<T,2>& data,const TV& A,const TV& B);

// Compute V(data_m,data_n); add to data
template<class T,class TV,int m,int n> void Combine_Data(DATA<T,4>& data,const DATA<T,2>& V,const DATA<TV,m>& data_m,const DATA<TV,n>& data_n,const VECTOR<int,m>& index_m,const VECTOR<int,n>& index_n);

// Individual cases
template<class T,class TV> void Case_CCAA(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_CCAB(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_CCBB(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_BCAC(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_BCBC(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_ACAC(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);

template<class T,class TV> void Area_From_Segments(DATA<T,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
}
}
#endif

