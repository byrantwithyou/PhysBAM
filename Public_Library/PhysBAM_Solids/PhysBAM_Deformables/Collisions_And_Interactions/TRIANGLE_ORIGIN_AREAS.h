//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_ORIGIN_AREAS
//##################################################################### 
#ifndef __TRIANGLE_ORIGIN_AREAS__
#define __TRIANGLE_ORIGIN_AREAS__

#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
namespace TRIANGLE_ORIGIN_AREAS
{
enum POINT_CASE {inside, beyond, outside};
template<class TV> POINT_CASE Classify_Point(const TV& A,const TV& B,const TV& C,const TV& P);

template<class T,int n>
struct VOL_DATA
{
    T V;
    VECTOR<T,3> G[n];
    MATRIX<T,3> H[n][n];
};

template<class T>
struct PT_DATA
{
    int n;
    int index[5];
    VECTOR<T,3> V;
    MATRIX<T,3> G[5];
    MATRIX<T,3> H[3][5][5];
};

template<class T,int n> void Clear(VOL_DATA<T,n>& data);

template<class T,class TV> void Data_From_Dof(PT_DATA<T>& data,const TV& A);

// Intersect ABC with OP
template<class T,class TV> void Intersect_Triangle_Point(PT_DATA<T>& data,const TV& A,const TV& B,const TV& C,const TV& P);

// Intersect ABC with PQ
template<class T,class TV> void Intersect_Triangle_Segment(PT_DATA<T>& data,const TV& A,const TV& B,const TV& C,const TV& P,const TV& Q);

// Intersect OAB with PQ
template<class T,class TV> void Intersect_Segment_Segment(PT_DATA<T>& data,const TV& A,const TV& B,const TV& P,const TV& Q);

// Volume of OABC
template<class T,class TV> void Volume_From_Points(PT_DATA<T>& data,const TV& A,const TV& B,const TV& C);

// Compute V(data_m,data_n); add to data
template<class T> void Combine_Data(VOL_DATA<T,6>& data,const VOL_DATA<T,3>& V,const PT_DATA<T>& data_m,const PT_DATA<T>& data_n,const PT_DATA<T>& data_p);

// Individual cases
// TODO: List the cases
//template<class T,class TV> void Case_CCCAAA(VOL_DATA<T>,const TV& A,const TV& B,const TV& C,const TV& D,const TV& E,const TV& F);
//template<class T,class TV> void Case_CCCBBB(VOL_DATA<T>,const TV& A,const TV& B,const TV& C,const TV& D,const TV& E,const TV& F);

template<class T,class TV> void Volume_From_Triangles(VOL_DATA<T,6>& data,TV A,TV B,TV C,TV D,TV E,TV F);
template<class T,class TV> void Volume_From_Triangles_Cut(VOL_DATA<T,6>& data,TV A,TV B,TV C,TV D,TV E,TV F);
}
}
#endif

