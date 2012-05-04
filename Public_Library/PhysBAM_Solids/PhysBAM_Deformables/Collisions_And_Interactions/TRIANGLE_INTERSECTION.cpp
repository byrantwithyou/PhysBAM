//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SEGMENT_ORIGIN_AREAS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_INTERSECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Triangle_Intersection_Area
//#####################################################################
template<class T,class TV> T PhysBAM::
Triangle_Intersection_Area(const TRIANGLE_2D<T>& a,const TRIANGLE_2D<T>& b,VECTOR<TV,6>& G,VECTOR<VECTOR<MATRIX<T,2>,6>,6>& H)
{
    T A=0;
    for(int i=0;i<3;i++){
        int in=(i+1)%3;
        for(int j=0;j<3;j++){
            int jn=(j+1)%3;
            VECTOR<int,4> I(i,in,j+3,jn+3);
            ORIGIN_AREAS::VOL_DATA<T,2,4> data;
            TV const X[]={a.X(i),a.X(in),b.X(j),b.X(jn)};
            ORIGIN_AREAS::Volume_From_Simplices(data,TV(),X);
            A+=data.V;

            for(int k=0;k<4;k++) G(I(k))+=(const TV&)data.G[k];
            for(int k=0;k<4;k++) for(int m=0;m<4;m++) H(I(k))(I(m))+=data.H[k][m];}}
    return A;
}
//#####################################################################
// Function Topology_Aware_Intersection_Test
//#####################################################################
template<class TV> bool PhysBAM::
Topology_Aware_Intersection_Test(VECTOR<int,3> a,VECTOR<int,3> b,ARRAY_VIEW<const TV> X)
{
    typedef typename TV::SCALAR T;
    bool b1=a.Contains(b(0)),b2=a.Contains(b(1)),b3=a.Contains(b(2));
    int c=b1+b2+b3;
    if(c==0) return SIMPLEX_INTERACTIONS<T>::Intersection(VECTOR<TV,3>(X.Subset(a)),VECTOR<TV,3>(X.Subset(b)));
    if(c>2) return false;

    T a1=TRIANGLE_2D<T>::Signed_Area(X(a.x),X(a.y),X(a.z)),a2=TRIANGLE_2D<T>::Signed_Area(X(b.x),X(b.y),X(b.z));
    if(!a1 || !a2) return false;
    if(c==2) return (a1>0)!=(a2>0);

    if(b3) cyclic_shift(b);
    else if(b2){cyclic_shift(b);cyclic_shift(b);}
    if(b.x!=a.x){cyclic_shift(a);if(b.x!=a.x){cyclic_shift(a);}}
    if(a1<0) exchange(a.y,a.z);
    if(a2<0) exchange(b.y,b.z);
    TV C(X(a.x));
    MATRIX<T,2> M1(X(a.y)-C,X(a.z)-C),M2(X(b.y)-C,X(b.z)-C);
    MATRIX<T,2> R=M1.Inverse()*M2;
    if(R.Column(0).All_Greater(TV()) || R.Column(1).All_Greater(TV())) return true;
    MATRIX<T,2> S=M2.Inverse()*M1;
    if(S.Column(0).All_Greater(TV()) || S.Column(1).All_Greater(TV())) return true;
	return false;
}
//#####################################################################
// Function Topology_Aware_Intersection_Test
//#####################################################################
template<class TV> bool PhysBAM::
Topology_Aware_Intersection_Test(VECTOR<int,4> a,VECTOR<int,4> b,ARRAY_VIEW<const TV> X)
{
    typedef typename TV::SCALAR T;
    bool b1=a.Contains(b(0)),b2=a.Contains(b(1)),b3=a.Contains(b(2)),b4=a.Contains(b(3));
    int c=b1+b2+b3+b4;
    if(c==0) return SIMPLEX_INTERACTIONS<T>::Intersection(VECTOR<TV,4>(X.Subset(a)),VECTOR<TV,4>(X.Subset(b)));
    if(c>3) return false;

    T a1=TETRAHEDRON<T>::Signed_Volume(X(a(0)),X(a(1)),X(a(2)),X(a(3))),a2=TETRAHEDRON<T>::Signed_Volume(X(b(0)),X(b(1)),X(b(2)),X(b(3)));
    if(!a1 || !a2) return false;
    if(c==3) return (a1>0)!=(a2>0);

    return true;
}
template bool PhysBAM::Topology_Aware_Intersection_Test<VECTOR<float,2> >(VECTOR<int,3>,VECTOR<int,3>,ARRAY_VIEW<VECTOR<float,2> const,int>);
template bool PhysBAM::Topology_Aware_Intersection_Test<VECTOR<float,3> >(VECTOR<int,4>,VECTOR<int,4>,ARRAY_VIEW<VECTOR<float,3> const,int>);
template float PhysBAM::Triangle_Intersection_Area<float,VECTOR<float,2> >(TRIANGLE_2D<float> const&,TRIANGLE_2D<float> const&,VECTOR<VECTOR<float,2>,6>&,VECTOR<VECTOR<MATRIX<float,2,2>,6>,6>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template bool PhysBAM::Topology_Aware_Intersection_Test<VECTOR<double,2> >(VECTOR<int,3>,VECTOR<int,3>,ARRAY_VIEW<VECTOR<double,2> const,int>);
template bool PhysBAM::Topology_Aware_Intersection_Test<VECTOR<double,3> >(VECTOR<int,4>,VECTOR<int,4>,ARRAY_VIEW<VECTOR<double,3> const,int>);
template double PhysBAM::Triangle_Intersection_Area<double,VECTOR<double,2> >(TRIANGLE_2D<double> const&,TRIANGLE_2D<double> const&,VECTOR<VECTOR<double,2>,6>&,VECTOR<VECTOR<MATRIX<double,2,2>,6>,6>&);
#endif
