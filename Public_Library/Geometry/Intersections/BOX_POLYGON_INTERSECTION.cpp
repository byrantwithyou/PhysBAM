//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/RANGE.h>
namespace PhysBAM{
namespace INTERSECTION{
template<class TV> void
Box_Polygon_Intersection(const RANGE<TV>& box,ARRAY<TV>& vertices,bool repeat_first_vertex)
{
    typedef typename TV::SCALAR T;
    if(!repeat_first_vertex) vertices.Append(TV(vertices(0)));

    TV bound[2]={box.min_corner,box.max_corner};
    ARRAY<TV> array;
    for(int s=0;s<2;s++){
        int sign=2*s-1;
        for(int a=0;a<2;a++){
            int first[2]={-1,-1};
            T b=bound[s](a);
            bool prev=(vertices(0)(a)-b)*sign>0; // 0=inside, 1=outside
            for(int i=1;i<vertices.m;i++){
                bool next=(vertices(i)(a)-b)*sign>0;
                if(next!=prev){
                    prev=next;
                    first[next]=i;}}
            
            // All in or all out
            if(first[0]<0){
                assert(first[1]<0);
                if(!prev) continue; // All in
                vertices.Remove_All(); // All out
                return;}
            assert(first[1]>=0);

            // P=enter, Q=leave
            TV A=vertices(first[0]-1),B=vertices(first[0]),P=bound[s];
            TV C=vertices(first[1]-1),D=vertices(first[1]),Q=bound[s];
            T iw=A(a)-b,ix=B(a)-b,jw=C(a)-b,jx=D(a)-b;
            T u=ix/(ix-iw),v=jx/(jx-jw);
            P(1-a)=u*A(1-a)+(1-u)*B(1-a);
            Q(1-a)=v*C(1-a)+(1-v)*D(1-a);

            if(first[0]<first[1]){ // ooiiioo -> PiiiQP
                array.Append(P);
                for(int i=first[0];i<first[1];i++)
                    array.Append(vertices(i));
                array.Append(Q);
                array.Append(P);}
            else{// iioooii -> iiQPii
                for(int i=0;i<first[1];i++)
                    array.Append(vertices(i));
                array.Append(Q);
                array.Append(P);
                for(int i=first[0];i<vertices.m;i++)
                    array.Append(vertices(i));}
            array.Exchange(vertices);
            array.Remove_All();}}

    if(!repeat_first_vertex && vertices.m) vertices.Pop();
}
template void Box_Polygon_Intersection<VECTOR<float,2> >(const RANGE<VECTOR<float,2> >& box,ARRAY<VECTOR<float,2> >& vertices,bool repeat_first_vertex);
template void Box_Polygon_Intersection<VECTOR<double,2> >(const RANGE<VECTOR<double,2> >& box,ARRAY<VECTOR<double,2> >& vertices,bool repeat_first_vertex);
}
}

