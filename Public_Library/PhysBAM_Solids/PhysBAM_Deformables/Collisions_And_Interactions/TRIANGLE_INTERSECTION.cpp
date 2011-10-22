//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_INTERSECTION.h>
using namespace PhysBAM;
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
            for(int k=1;k<=4;k++) for(int m=1;m<=4;m++) G(I(k),I(m))+=tG(k,m);
        }
    }
}


