//#####################################################################
// Copyright 2007, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXACT_SIMPLEX_INTERACTIONS
//##################################################################### 
#include <Deformables/Collisions_And_Interactions/EXACT_SIMPLEX_INTERACTIONS.h>
using namespace PhysBAM;
//#####################################################################
// Function Exact_Segment_Intersection
//#####################################################################
bool EXACT_SIMPLEX_INTERACTIONS::
Exact_Segment_Intersection(const VECTOR<float,2>& p1,const VECTOR<float,2>& p2,const VECTOR<float,2>& q3,const VECTOR<float,2>& q4,
    const int rank1,const int rank2,const int rank3,const int rank4)
{
    const VECTOR<VECTOR<float,2>,4> points(p1,p2,q3,q4);
    const VECTOR<int,4> ranks(rank1,rank2,rank3,rank4);
    bool psa[16];
    int plus,minus;
    unsigned int mask;

    // Compute disambiguated signed areas for all triangles formed by three input points
    for(unsigned int i=0;i<2;i++) for(unsigned int j=i+1;j<2;j++) for(unsigned int k=j+1;k<3;k++)
        psa[(1<<i)|(1<<j)|(1<<k)]=Positive_Signed_Area(points(i+1),points(j+1),points(k+1),ranks(i+1),ranks(j+1),ranks(k+1));

    // Separator line must contain a vertex from each segment
    for(unsigned int i=0;i<2;i++) for(unsigned int j=1;j<3;j++){
        plus=minus=0;mask=(1<<i)|(1<<j);
        for(unsigned int k=0;k<i;k++)    if(psa[(1<<k)|mask]) plus++;  else minus++;
        for(unsigned int k=i+1;k<2;k++) if(psa[(1<<k)|mask]) minus++; else plus++;
        for(unsigned int k=2;k<j;k++)    if(psa[(1<<k)|mask]) plus++;  else minus++;
        for(unsigned int k=j+1;k<4;k++) if(psa[(1<<k)|mask]) minus++; else plus++;
        if(plus==0||minus==0) return false;}
    // No separator line found, objects are intersecting
    return true;
}
//#####################################################################
// Function Exact_Segment_Intersection
//#####################################################################
// Uses lexicographical ordering of vertices (assumes no vertices are coincident)
bool EXACT_SIMPLEX_INTERACTIONS::
Exact_Segment_Intersection(const VECTOR<float,2>& p1,const VECTOR<float,2>& p2,const VECTOR<float,2>& q3,const VECTOR<float,2>& q4)
{
    const VECTOR<VECTOR<float,2>,4> points(p1,p2,q3,q4);
    bool psa[16];
    int plus,minus;
    unsigned int mask;

    // Compute disambiguated signed areas for all triangles formed by three input points
    for(unsigned int i=0;i<2;i++) for(unsigned int j=i+1;j<2;j++) for(unsigned int k=j+1;k<3;k++) psa[(1<<i)|(1<<j)|(1<<k)]=Positive_Signed_Area(points(i+1),points(j+1),points(k+1));

    // Separator line must contain a vertex from each segment
    for(unsigned int i=0;i<2;i++) for(unsigned int j=1;j<3;j++){
        plus=minus=0;mask=(1<<i)|(1<<j);
        for(unsigned int k=0;k<i;k++)    if(psa[(1<<k)|mask]) plus++;  else minus++;
        for(unsigned int k=i+1;k<2;k++) if(psa[(1<<k)|mask]) minus++; else plus++;
        for(unsigned int k=2;k<j;k++)    if(psa[(1<<k)|mask]) plus++;  else minus++;
        for(unsigned int k=j+1;k<4;k++) if(psa[(1<<k)|mask]) minus++; else plus++;
        if(plus==0||minus==0) return false;}
    // No separator line found, objects are intersecting
    return true;
}
//#####################################################################
// Function Positive_Signed_Area
//#####################################################################
// Using Simulation Of Simplicity
bool EXACT_SIMPLEX_INTERACTIONS::
Positive_Signed_Area(VECTOR<float,2> x0,VECTOR<float,2> x1,VECTOR<float,2> x2,int rank1,int rank2,int rank3)
{
    bool positive=true;
    if(rank1>rank2){exchange(rank1,rank2);exchange(x0,x1);positive=!positive;}
    if(rank1>rank3){exchange(rank1,rank3);exchange(x0,x2);positive=!positive;}
    if(rank2>rank3){exchange(rank2,rank3);exchange(x1,x2);positive=!positive;}
    double area=Exact_Signed_Area(x0,x1,x2);
    if(area!=0) return (area<0)^positive; // works with area=-0
    else if(x2.x!=x1.x) return (x2.x<x1.x)^positive;
    else if(x1.y!=x2.y) return (x1.y<x2.y)^positive;
    else if(x0.x!=x2.x) return (x0.x<x2.x)^positive;
    return positive;
}
//#####################################################################
// Function Positive_Signed_Area
//#####################################################################
// Using Simulation Of Simplicity and lexicographical ordering of vertices (assumes no vertices are coincident)
bool EXACT_SIMPLEX_INTERACTIONS::
Positive_Signed_Area(VECTOR<float,2> x0,VECTOR<float,2> x1,VECTOR<float,2> x2)
{
    bool positive=true;
    if(x0.x>x1.x || (x0.x==x1.x && x0.y>x1.y)){exchange(x0,x1);positive=!positive;}
    if(x0.x>x2.x || (x0.x==x2.x && x0.y>x2.y)){exchange(x0,x2);positive=!positive;}
    if(x1.x>x2.x || (x1.x==x2.x && x1.y>x2.y)){exchange(x1,x2);positive=!positive;}
    double area=Exact_Signed_Area(x0,x1,x2);
    if(area!=0) return (area<0)^positive; // works with area=-0
    if(x2.x!=x1.x) return !positive;
    return positive;
}
//#####################################################################
// Function Exact_Signed_Area
//#####################################################################
// Sign is guaranteed to be exact (if area is zero, we get exactly zero)
double EXACT_SIMPLEX_INTERACTIONS::
Exact_Signed_Area(const VECTOR<float,2>& x0,const VECTOR<float,2>& x1,const VECTOR<float,2>& x2)
{
    VECTOR<double,6> terms;
    terms(0)=((double)x0.x)*((double)x1.y);terms(1)=((double)x1.x)*((double)x2.y);terms(2)=((double)x2.x)*((double)x0.y);
    terms(3)=-((double)x0.y)*((double)x1.x);terms(4)=-((double)x1.y)*((double)x2.x);terms(5)=-((double)x2.y)*((double)x0.x);
    for(int i=0;i<5;i++) for(int j=i+1;j<6;j++) if(abs(terms(i))<abs(terms(j))) exchange(terms(i),terms(j));
    for(int number_of_terms=4;number_of_terms>=0;number_of_terms--){
        int i=1;for(;i<number_of_terms-1;i++) if(terms(0)*terms(i)<0) break;
        terms(i-1)+=terms(i);for(int j=i;j<number_of_terms-1;j++) terms(j)=terms(j+1);
        for(int j=i-1;j<number_of_terms-2;j++) if(abs(terms(j))<abs(terms(j+1))) exchange(terms(j),terms(j+1));}
    return terms(0);
}
//####################################################################

