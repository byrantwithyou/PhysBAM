//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRAPEZOID_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_INTERSECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Triangle_Intersection_Area
//#####################################################################
template<class T,class TV> T PhysBAM::
Triangle_Intersection_Area(const TRIANGLE_2D<T>& a,const TRIANGLE_2D<T>& b,VECTOR<TV,6>& G,VECTOR<VECTOR<MATRIX<T,2>,6>,6>& H)
{
    T A=0;
    for(int i=1;i<=3;i++){
        int in=i%3+1;
        for(int j=1;j<=3;j++){
            int jn=j%3+1;
            VECTOR<TV,4> tG;
            VECTOR<VECTOR<MATRIX<T,2>,4>,4> tH;
            A+=Trapezoid_Intersection_Area(a.X(i),a.X(in),b.X(j),b.X(jn),tG,tH);
            VECTOR<int,4> I(i,in,j,in);
            for(int k=1;k<=4;k++) G(I(k))=tG(k);
            for(int k=1;k<=4;k++) for(int m=1;m<=4;m++) G(I(k))(I(m))+=tG(k)(m);
        }
    }
    return A;
}
//#####################################################################
// Function Topology_Aware_Triangle_Intersection_Test
//#####################################################################
template<class TV> bool PhysBAM::
Topology_Aware_Triangle_Intersection_Test(VECTOR<int,3> a,VECTOR<int,3> b,ARRAY_VIEW<const TV> X)
{
    typedef typename TV::SCALAR T;
    bool b1=a.Contains(b(1)),b2=a.Contains(b(2)),b3=a.Contains(b(3));
    int c=b1+b2+b3;
    if(c==0) return SIMPLEX_INTERACTIONS<T>::Intersection(VECTOR<TV,3>(X.Subset(a)),VECTOR<TV,3>(X.Subset(b)));
    if(c>2) return false;

    T a1=TRIANGLE_2D<T>::Signed_Area(X(a.x),X(a.y),X(a.z)),a2=TRIANGLE_2D<T>::Signed_Area(X(b.x),X(b.y),X(b.z));
    if(!a1 || !a2) return false;
    if(c==2) return (a1>0)!=(a2>0);

    if(a1<0) exchange(a.x,a.y);
    if(a2<0) exchange(b.x,b.y);
    if(b2) cyclic_shift(a);
    else if(b3){cyclic_shift(a);cyclic_shift(a);}
    if(b.x!=a.x){cyclic_shift(b);if(b.x!=a.x){cyclic_shift(b);}}
    TV C(X(a.x));
    MATRIX<T,2> M1(X(a.y)-C,X(a.z)-C),M2(X(b.y)-C,X(b.z)-C);
    MATRIX<T,2> R=M1.Cofactor_Matrix()*M2;
    return R.Column(1).All_Greater(TV()) || R.Column(2).All_Greater(TV());
}
template bool Topology_Aware_Triangle_Intersection_Test<VECTOR<float,2> >(VECTOR<int,3>,VECTOR<int,3>,ARRAY_VIEW<VECTOR<float,2> const,int>);
template float Triangle_Intersection_Area<float,VECTOR<float,2> >(TRIANGLE_2D<float> const&,TRIANGLE_2D<float> const&,VECTOR<VECTOR<float,2>,6>&,VECTOR<VECTOR<MATRIX<float,2,2>,6>,6>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool Topology_Aware_Triangle_Intersection_Test<VECTOR<double,2> >(VECTOR<int,3>,VECTOR<int,3>,ARRAY_VIEW<VECTOR<double,2> const,int>);
template double Triangle_Intersection_Area<double,VECTOR<double,2> >(TRIANGLE_2D<double> const&,TRIANGLE_2D<double> const&,VECTOR<VECTOR<double,2>,6>&,VECTOR<VECTOR<MATRIX<double,2,2>,6>,6>&);
#endif
