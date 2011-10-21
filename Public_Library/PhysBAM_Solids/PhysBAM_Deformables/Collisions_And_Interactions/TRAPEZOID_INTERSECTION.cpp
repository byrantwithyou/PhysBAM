//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRAPEZOID_INTERSECTION.h>
using namespace PhysBAM;
// Case 1 order: a c b d; ou = c over, b under
template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1ou(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    LOG::cout<<__FUNCTION__<<std::endl;
    return (T).5*(-d.x+c.x)*(-2*b.y*sqr(c.x)*a.y-2*d.y*d.y*b.x*a.x+b.y*b.y*c.x*c.x+a.y*a.y*c.x*c.x+d.y*d.y*b.x*b.x+d.y*d.y*a.x*a.x-b.y*d.x*d.y*b.x+b.y*d.x*d.y*a.x-b.y*d.x*c.y*b.x+b.y*d.x*c.y*a.x+b.y*c.x*d.y*b.x-b.y*c.x*d.y*a.x-b.y*c.x*c.y*b.x+b.y*c.x*c.y*a.x+a.y*d.x*d.y*b.x-a.y*d.x*d.y*a.x+a.y*d.x*c.y*b.x-a.y*d.x*c.y*a.x-a.y*c.x*d.y*b.x+a.y*c.x*d.y*a.x+a.y*c.x*c.y*b.x-a.y*c.x*c.y*a.x-2*a.x*b.y*b.y*c.x-2*a.y*c.y*b.x*b.x+2*b.y*b.x*c.y*a.x+2*a.y*b.x*b.y*c.x+2*a.y*b.x*c.y*a.x+2*a.y*a.x*b.y*c.x-2*a.y*a.y*c.x*b.x-2*a.y*b.x*a.x*b.y-2*b.y*a.x*a.x*c.y+a.y*a.y*b.x*b.x+a.x*a.x*b.y*b.y)/(-b.x+a.x)/(-b.y*d.x+b.y*c.x+a.y*d.x-a.y*c.x+d.y*b.x-d.y*a.x-c.y*b.x+c.y*a.x);
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1oo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    LOG::cout<<__FUNCTION__<<std::endl;
    return -(T).5*(-c.y*b.y*c.x*c.x-d.y*b.y*c.x*c.x+c.y*a.y*c.x*c.x+d.y*a.y*c.x*c.x+b.x*b.y*b.y*c.x-2*a.y*c.x*d.y*b.x+c.y*c.y*c.x*b.x-2*c.y*c.x*a.y*d.x+2*c.y*c.x*b.y*d.x+2*b.y*c.x*d.y*a.x-a.x*b.y*b.y*c.x-c.y*c.y*c.x*a.x-a.y*c.y*b.x*b.x+a.y*d.y*b.x*b.x+2*b.y*b.x*c.y*a.x+b.y*d.y*b.x*b.x+a.x*b.y*b.y*d.x-2*b.y*b.x*d.y*a.x-2*b.y*d.x*c.y*a.x+2*a.y*d.x*c.y*b.x-b.x*b.y*b.y*d.x-c.y*c.y*b.x*d.x-b.y*c.y*b.x*b.x+c.y*c.y*a.x*d.x)/(-b.y*d.x+b.y*c.x+a.y*d.x-a.y*c.x+d.y*b.x-d.y*a.x-c.y*b.x+c.y*a.x);
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1uu(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    LOG::cout<<__FUNCTION__<<std::endl;
    return -(T).5*(4*a.y*b.x*d.x*a.x*b.y*c.x+2*a.y*b.x*d.x*c.x*d.y*a.x-2*a.y*c.x*b.x*c.y*a.x*d.x-2*a.x*b.y*d.x*c.x*d.y*b.x-6*a.x*b.y*c.x*c.y*b.x*d.x+4*c.y*b.x*d.x*c.x*d.y*a.x+4*a.x*b.y*d.x*d.x*c.y*b.x+4*a.x*a.x*b.y*c.x*c.y*d.x+2*a.x*b.y*c.x*c.x*d.y*b.x-2*a.y*b.x*d.x*d.x*a.x*b.y-2*a.y*b.x*b.x*d.x*c.x*d.y-2*a.y*c.x*c.x*b.x*a.x*b.y+2*a.y*c.x*b.x*b.x*c.y*d.x-2*a.y*c.x*c.x*b.x*d.y*a.x+a.y*a.y*b.x*b.x*d.x*d.x+a.y*a.y*c.x*c.x*b.x*b.x+a.x*a.x*b.y*b.y*d.x*d.x+a.x*a.x*b.y*b.y*c.x*c.x+c.y*c.y*b.x*b.x*d.x*d.x+c.y*c.y*a.x*a.x*d.x*d.x+c.x*c.x*d.y*d.y*b.x*b.x+c.x*c.x*d.y*d.y*a.x*a.x-2*b.x*b.x*b.x*d.y*d.y*a.x-2*b.x*b.x*b.x*b.x*d.y*c.y-2*b.x*b.x*b.x*c.y*c.y*a.x+b.x*b.x*d.y*d.y*a.x*a.x+b.x*b.x*c.y*c.y*a.x*a.x+b.x*b.x*b.x*b.x*d.y*d.y+b.x*b.x*b.x*b.x*c.y*c.y+c.x*c.x*c.x*c.x*b.y*b.y+c.x*c.x*c.x*c.x*a.y*a.y-2*a.y*a.y*b.x*b.x*d.x*c.x+2*a.y*c.x*c.x*b.x*b.x*d.y-2*a.x*a.x*b.y*b.y*d.x*c.x-2*a.x*a.x*b.y*d.x*d.x*c.y-2*c.y*c.y*b.x*d.x*d.x*a.x-2*c.x*c.x*d.y*d.y*b.x*a.x-2*b.x*b.x*d.y*a.x*a.x*c.y-b.x*b.x*b.x*b.y*d.x*d.y+b.x*b.x*b.x*b.y*d.x*c.y+b.x*b.x*b.x*b.y*c.x*d.y-b.x*b.x*b.x*b.y*c.x*c.y+b.x*b.x*b.x*a.y*d.x*d.y-b.x*b.x*b.x*a.y*d.x*c.y-b.x*b.x*b.x*a.y*c.x*d.y+b.x*b.x*b.x*a.y*c.x*c.y+4*b.x*b.x*b.x*d.y*c.y*a.x-2*c.x*c.x*b.y*d.x*d.x*a.y+4*c.x*c.x*c.x*b.y*d.x*a.y+c.x*c.x*c.x*b.y*d.y*b.x-c.x*c.x*c.x*b.y*d.y*a.x-c.x*c.x*c.x*b.y*c.y*b.x+c.x*c.x*c.x*b.y*c.y*a.x-c.x*c.x*c.x*a.y*d.y*b.x+c.x*c.x*c.x*a.y*d.y*a.x+c.x*c.x*c.x*a.y*c.y*b.x-c.x*c.x*c.x*a.y*c.y*a.x-2*c.y*b.x*b.x*d.x*d.x*b.y+2*c.y*b.x*b.x*b.x*d.x*d.y+4*c.y*c.y*b.x*b.x*d.x*a.x-2*c.y*c.y*b.x*d.x*a.x*a.x-2*c.x*c.x*d.y*b.x*b.x*b.y+4*c.x*d.y*d.y*b.x*b.x*a.x+2*c.x*d.y*b.x*b.x*b.x*c.y-2*c.x*d.y*d.y*b.x*a.x*a.x-2*c.y*b.x*b.x*d.x*c.x*d.y-2*c.y*a.x*a.x*d.x*c.x*d.y+b.x*b.x*b.y*d.x*d.y*a.x-b.x*b.x*b.y*d.x*c.y*a.x-b.x*b.x*b.y*c.x*d.y*a.x+b.x*b.x*b.y*c.x*c.y*a.x-b.x*b.x*a.y*d.x*d.y*a.x+b.x*b.x*a.y*d.x*c.y*a.x+b.x*b.x*a.y*c.x*d.y*a.x-b.x*b.x*a.y*c.x*c.y*a.x-c.x*c.x*b.y*d.x*d.y*b.x+c.x*c.x*b.y*d.x*d.y*a.x+c.x*c.x*b.y*d.x*c.y*b.x-c.x*c.x*b.y*d.x*c.y*a.x+c.x*c.x*a.y*d.x*d.y*b.x-c.x*c.x*a.y*d.x*d.y*a.x-c.x*c.x*a.y*d.x*c.y*b.x+c.x*c.x*a.y*d.x*c.y*a.x+2*c.y*b.x*b.x*d.x*b.y*c.x-4*c.y*b.x*b.x*d.x*d.y*a.x+2*c.y*b.x*d.x*a.x*a.x*d.y+2*c.x*d.y*b.x*b.x*b.y*d.x-4*c.x*d.y*b.x*b.x*c.y*a.x+2*c.x*d.y*b.x*a.x*a.x*c.y-2*c.x*c.x*c.x*b.y*b.y*d.x-2*c.x*c.x*c.x*c.x*b.y*a.y-2*c.x*c.x*c.x*a.y*a.y*d.x+c.x*c.x*b.y*b.y*d.x*d.x+c.x*c.x*a.y*a.y*d.x*d.x-2*c.y*c.y*b.x*b.x*b.x*d.x-2*c.x*d.y*d.y*b.x*b.x*b.x-2*b.y*b.y*c.x*c.x*c.x*a.x-2*a.y*a.y*c.x*c.x*c.x*b.x-2*c.x*b.y*b.y*d.x*d.x*a.x+4*c.x*c.x*b.y*b.y*d.x*a.x+2*b.y*c.x*c.x*c.x*b.x*a.y+2*b.y*c.x*c.x*c.x*a.x*a.y-2*b.y*c.x*c.x*a.x*a.x*c.y-2*c.x*a.y*a.y*d.x*d.x*b.x+4*c.x*c.x*a.y*a.y*d.x*b.x-2*a.y*c.x*c.x*b.x*b.x*c.y+2*c.x*b.y*d.x*d.x*b.x*a.y-4*c.x*c.x*b.y*d.x*b.x*a.y+2*c.x*b.y*d.x*d.x*a.x*a.y-4*c.x*c.x*b.y*d.x*a.x*a.y+2*b.y*c.x*c.x*b.x*c.y*a.x+2*a.y*c.x*c.x*b.x*c.y*a.x)/(-d.x+c.x)/(-b.x+a.x)/(-b.y*d.x+b.y*c.x+a.y*d.x-a.y*c.x+d.y*b.x-d.y*a.x-c.y*b.x+c.y*a.x);
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1uo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    LOG::cout<<__FUNCTION__<<std::endl;
    return -(T).5*(-d.x+c.x)*(-2*b.y*d.x*d.x*a.y-2*c.y*c.y*b.x*a.x+b.y*b.y*d.x*d.x+a.y*a.y*d.x*d.x+c.y*c.y*b.x*b.x+c.y*c.y*a.x*a.x-b.y*d.x*d.y*b.x+b.y*d.x*d.y*a.x+b.y*d.x*c.y*b.x-b.y*d.x*c.y*a.x-b.y*c.x*d.y*b.x+b.y*c.x*d.y*a.x-b.y*c.x*c.y*b.x+b.y*c.x*c.y*a.x+a.y*d.x*d.y*b.x-a.y*d.x*d.y*a.x-a.y*d.x*c.y*b.x+a.y*d.x*c.y*a.x+a.y*c.x*d.y*b.x-a.y*c.x*d.y*a.x+a.y*c.x*c.y*b.x-a.y*c.x*c.y*a.x-2*a.x*b.y*b.y*d.x-2*a.y*d.y*b.x*b.x+2*b.y*b.x*d.y*a.x+2*a.y*b.x*b.y*d.x+2*a.y*b.x*d.y*a.x+2*a.y*a.x*b.y*d.x-2*a.y*a.y*b.x*d.x-2*a.y*b.x*a.x*b.y+a.y*a.y*b.x*b.x+a.x*a.x*b.y*b.y-2*a.x*a.x*b.y*d.y)/(-b.x+a.x)/(-b.y*d.x+b.y*c.x+a.y*d.x-a.y*c.x+d.y*b.x-d.y*a.x-c.y*b.x+c.y*a.x);
}

// Case 2 order: a c d b; ou = c over, d under
template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2ou(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    LOG::cout<<__FUNCTION__<<std::endl;
    return (T).5*(-d.x+c.x)*(-2*b.y*c.x*c.x*a.y-2*d.y*d.y*b.x*a.x+b.y*b.y*c.x*c.x+a.y*a.y*c.x*c.x+d.y*d.y*b.x*b.x+d.y*d.y*a.x*a.x-b.y*d.x*d.y*b.x+b.y*d.x*d.y*a.x-b.y*d.x*c.y*b.x+b.y*d.x*c.y*a.x+b.y*c.x*d.y*b.x-b.y*c.x*d.y*a.x-b.y*c.x*c.y*b.x+b.y*c.x*c.y*a.x+a.y*d.x*d.y*b.x-a.y*d.x*d.y*a.x+a.y*d.x*c.y*b.x-a.y*d.x*c.y*a.x-a.y*c.x*d.y*b.x+a.y*c.x*d.y*a.x+a.y*c.x*c.y*b.x-a.y*c.x*c.y*a.x-2*a.x*b.y*b.y*c.x-2*a.y*c.y*b.x*b.x+2*b.y*b.x*c.y*a.x+2*a.y*b.x*b.y*c.x+2*a.y*b.x*c.y*a.x+2*a.y*a.x*b.y*c.x-2*a.y*a.y*c.x*b.x-2*a.y*b.x*a.x*b.y-2*b.y*a.x*a.x*c.y+a.y*a.y*b.x*b.x+a.x*a.x*b.y*b.y)/(-b.x+a.x)/(-b.y*d.x+b.y*c.x+a.y*d.x-a.y*c.x+d.y*b.x-d.y*a.x-c.y*b.x+c.y*a.x);
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2oo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    LOG::cout<<__FUNCTION__<<std::endl;
    return -(T).5*(-d.x+c.x)*(-b.y*c.x+a.y*c.x-2*a.y*b.x+2*a.x*b.y-b.y*d.x+a.y*d.x)/(-b.x+a.x);
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2uu(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    LOG::cout<<__FUNCTION__<<std::endl;
    return -(T).5*(d.y+c.y)*(-d.x+c.x);
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2uo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    LOG::cout<<__FUNCTION__<<std::endl;
    return -(T).5*(-d.x+c.x)*(-2*b.y*d.x*d.x*a.y-2*c.y*c.y*b.x*a.x+b.y*b.y*d.x*d.x+a.y*a.y*d.x*d.x+c.y*c.y*b.x*b.x+c.y*c.y*a.x*a.x-b.y*d.x*d.y*b.x+b.y*d.x*d.y*a.x+b.y*d.x*c.y*b.x-b.y*d.x*c.y*a.x-b.y*c.x*d.y*b.x+b.y*c.x*d.y*a.x-b.y*c.x*c.y*b.x+b.y*c.x*c.y*a.x+a.y*d.x*d.y*b.x-a.y*d.x*d.y*a.x-a.y*d.x*c.y*b.x+a.y*d.x*c.y*a.x+a.y*c.x*d.y*b.x-a.y*c.x*d.y*a.x+a.y*c.x*c.y*b.x-a.y*c.x*c.y*a.x-2*a.x*b.y*b.y*d.x-2*a.y*d.y*b.x*b.x+2*b.y*b.x*d.y*a.x+2*a.y*b.x*b.y*d.x+2*a.y*b.x*d.y*a.x+2*a.y*a.x*b.y*d.x-2*a.y*a.y*b.x*d.x-2*a.y*b.x*a.x*b.y+a.y*a.y*b.x*b.x+a.x*a.x*b.y*b.y-2*a.x*a.x*b.y*d.y)/(-b.x+a.x)/(-b.y*d.x+b.y*c.x+a.y*d.x-a.y*c.x+d.y*b.x-d.y*a.x-c.y*b.x+c.y*a.x);
}

// Case 1 order: a c b d; ou = c over, b under
template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    bool co=(-a.y*b.x-b.y*c.x+a.x*b.y+a.y*c.x+c.y*b.x-c.y*a.x)>0;
    bool bo=(b.y*d.x-b.y*c.x-c.y*d.x-d.y*b.x+c.x*d.y+c.y*b.x)>0;
    LOG::cout<<co<<bo<<std::endl;
    if(co){
        if(bo) return Trapezoid_Intersection_Area_Case_1oo(a,b,c,d,dA,H);
        else return Trapezoid_Intersection_Area_Case_1ou(a,b,c,d,dA,H);}
    else{
        if(bo) return Trapezoid_Intersection_Area_Case_1uo(a,b,c,d,dA,H);
        else return Trapezoid_Intersection_Area_Case_1uu(a,b,c,d,dA,H);}
}

// Case 2 order: a c d b; ou = c over, d under
template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    bool co=(-a.y*b.x-b.y*c.x+a.x*b.y+a.y*c.x+c.y*b.x-c.y*a.x)>0;
    bool Do=d.y*b.x-d.y*a.x-a.y*b.x-b.y*d.x+a.x*b.y+a.y*d.x;
    LOG::cout<<co<<Do<<std::endl;
    if(co){
        if(Do) return Trapezoid_Intersection_Area_Case_2oo(a,b,c,d,dA,H);
        else return Trapezoid_Intersection_Area_Case_2ou(a,b,c,d,dA,H);}
    else{
        if(Do) return Trapezoid_Intersection_Area_Case_2uo(a,b,c,d,dA,H);
        else return Trapezoid_Intersection_Area_Case_2uu(a,b,c,d,dA,H);}
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area(TV a,TV b,TV c,TV d,VECTOR<TV,4>& dA,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    VECTOR<int,4> indices(1,2,3,4);
    T sign=1,A=0;
    if(a.x>b.x){sign=-sign;exchange(a.x,b.x);exchange(indices(1),indices(2));}
    if(c.x>d.x){sign=-sign;exchange(c.x,d.x);exchange(indices(3),indices(4));}
    if(a.x>c.x){exchange(a.x,c.x);exchange(b.x,d.x);exchange(indices(1),indices(3));exchange(indices(2),indices(4));}
    LOG::cout<<indices<<std::endl;

    VECTOR<TV,4> tdA;
    VECTOR<VECTOR<MATRIX<T,2>,4>,4> tH;

    if(b.x<c.x){dA=tdA;H=tH;return 0;}

    
    if(b.x<d.x) A=sign*Trapezoid_Intersection_Area_Case_1(a,b,c,d,tdA,tH); // a c b d
    else A=sign*Trapezoid_Intersection_Area_Case_2(a,b,c,d,tdA,tH); // a c d b
    for(int i=1;i<=4;i++){
        dA(indices(i))=sign*tdA(i);
        for(int j=1;j<=4;j++) H(indices(i))(indices(j))=sign*tH(i)(j);}
    return A;
}

template float Trapezoid_Intersection_Area<float,VECTOR<float,2> >(VECTOR<float,2>,VECTOR<float,2>,VECTOR<float,2>,VECTOR<float,2>,VECTOR<VECTOR<float,2>,4>&,VECTOR<VECTOR<MATRIX<float,2,2>,4>,4>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template double Trapezoid_Intersection_Area<double,VECTOR<double,2> >(VECTOR<double,2>,VECTOR<double,2>,VECTOR<double,2>,VECTOR<double,2>,VECTOR<VECTOR<double,2>,4>&,VECTOR<VECTOR<MATRIX<double,2,2>,4>,4>&);
#endif
