//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_ORIGIN_AREAS
//##################################################################### 
#ifndef __SEGMENT_ORIGIN_AREAS__
#define __SEGMENT_ORIGIN_AREAS__

#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
namespace SEGMENT_ORIGIN_AREAS
{

template<class T,int m,int n>
struct DATA
{
    VECTOR<T,m> V;
    MATRIX<T,m,2> G[n];
    MATRIX<T,2> H[m][n][n];
};

template<class T,int m,int n> void Clear(DATA<T,m,n>& data);

template<class T,class TV> void Data_From_Dof(DATA<T,2,1>& data,const TV& A);

// Intersect AB with OP
template<class T,class TV> void Intersect_Segment_Point(DATA<T,2,3>& data,const TV& A,const TV& B,const TV& P);

// Intersect AB with CD
template<class T,class TV> void Intersect_Segments(DATA<T,2,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);

// Area of OAB
template<class T,class TV> void Area_From_Points(DATA<T,1,2>& data,const TV& A,const TV& B);

// Compute V(data_m,data_n); add to data
template<class T,int m,int n> void Combine_Data(DATA<T,1,4>& data,const DATA<T,1,2>& V,const DATA<T,2,m>& data_m,const DATA<T,2,n>& data_n,const int index_m[m],const int index_n[n]);

// Individual cases
template<class T,class TV> void Case_CCAA(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_CCAB(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_CCBB(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_BCAC(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_BCBC(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);
template<class T,class TV> void Case_ACAC(DATA<T,1,4>& data,const TV& A,const TV& B,const TV& C,const TV& D);

template<class T,class TV> void Area_From_Segments(DATA<T,1,4>& data,TV A,TV B,TV C,TV D);
}
}
#endif

