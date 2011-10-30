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

template<class T,int m,int n>
struct DATA
{
    VECTOR<T,m> V;
    MATRIX<T,m,3> G[n];
    MATRIX<T,3> H[m][n][n];
};

template<class T,int m,int n> void Clear(DATA<T,m,n>& data);

template<class T,class TV> void Data_From_Dof(DATA<T,3,1>& data,const TV& A);

// Intersect ABC with OP
template<class T,class TV> void Intersect_Triangle_Point(DATA<T,3,4>& data,const TV& A,const TV& B,const TV& C,const TV& P);

// Intersect ABC with PQ
template<class T,class TV> void Intersect_Triangle_Segment(DATA<T,3,5>& data,const TV& A,const TV& B,const TV& C,const TV& P,const TV& Q);

// Intersect OAB with PQ
template<class T,class TV> void Intersect_Segment_Segment(DATA<T,3,4>& data,const TV& A,const TV& B,const TV& P,const TV& Q);

// Area of OABC
template<class T,class TV> void Area_From_Points(DATA<T,1,3>& data,const TV& A,const TV& B,const TV& C);

// Compute V(data_m,data_n); add to data
template<class T,int m,int n,int p> void Combine_Data(DATA<T,1,6>& data,const DATA<T,1,3>& V,const DATA<T,3,m>& data_m,const DATA<T,3,n>& data_n,const DATA<T,3,n>& data_p,
    const int index_m[m],const int index_n[n],const int index_p[p]);

// Individual cases
// TODO: List the cases
template<class T,class TV> void Case_CCCAAA(DATA<T,1,6>& data,const TV& A,const TV& B,const TV& C,const TV& D,const TV& E,const TV& F);

template<class T,class TV> void Volume_From_Triangles(DATA<T,1,6>& data,TV A,TV B,TV C,TV D,TV E,TV F);
}
}
#endif

