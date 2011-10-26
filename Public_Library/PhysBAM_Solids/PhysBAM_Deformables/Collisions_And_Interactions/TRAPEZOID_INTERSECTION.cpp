//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRAPEZOID_INTERSECTION.h>
using namespace PhysBAM;
// Case 1 order: a c b d; ou = c over, b under
//extern ARRAY<int> trap_cases;
template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1ou(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
//    trap_cases.Append(1);
//    LOG::cout<<__FUNCTION__<<std::endl;
    T yvab=(a.y+b.y)/2,yba=b.y-a.y,yvcd=(c.y+d.y)/2,ydc=d.y-c.y,xba=b.x-a.x,xdc=d.x-c.x,xca=c.x-a.x;(void)yvab;(void)yba;(void)yvcd;(void)ydc;(void)xba;(void)xdc;(void)xca;
    T xbc=b.x-c.x,xbc_xba=xbc/xba,xca_xba=xca/xba,yba_xba=yba/xba,xca_xba2=xca_xba*xca_xba,xbc_xba2=xbc_xba*xbc_xba;

    T A=(T).5*xbc*(2*yvab+xca_xba*yba);
    G(1)(1)=-(T).5*xbc_xba2*yba;
    G(1)(2)=(T).5*xbc*xbc_xba;
    G(2)(1)=(T).5*(yba*xca_xba2+2*yvab);
    G(2)(2)=(T).5*(xba-xca*xca_xba);
    G(3)(1)=(T).5*(-2*yvab+yba-2*yba*xca_xba);
    G(3)(2)=0;
    G(4)(1)=0;
    G(4)(2)=0;
    H(1)(1)(1,1)=-xbc_xba2*yba_xba;
    H(1)(1)(2,1)=H(1)(1)(1,2)=(T).5*xbc_xba2;
    H(1)(1)(2,2)=0;
    H(2)(1)(1,1)=H(1)(2)(1,1)=-xbc_xba*xca_xba*yba_xba;
    H(2)(1)(1,2)=H(1)(2)(2,1)=(T).5*(1-xca_xba2);
    H(2)(2)(1,1)=-xca_xba2*yba_xba;
    H(2)(1)(2,1)=H(1)(2)(1,2)=-(T).5*xbc_xba2;
    H(2)(1)(2,2)=H(1)(2)(2,2)=0;
    H(2)(2)(2,1)=H(2)(2)(1,2)=(T).5*(1+xca_xba2);
    H(2)(2)(2,2)=0;
    H(3)(1)(1,1)=H(1)(3)(1,1)=xbc_xba*yba_xba;
    H(3)(1)(1,2)=H(1)(3)(2,1)=-xbc_xba;
    H(3)(2)(1,1)=H(2)(3)(1,1)=yba_xba*xca_xba;
    H(3)(2)(1,2)=H(2)(3)(2,1)=-xca_xba;
    H(3)(3)(1,1)=-yba_xba;
    H(3)(1)(2,1)=H(1)(3)(1,2)=0;
    H(3)(1)(2,2)=H(1)(3)(2,2)=0;
    H(3)(2)(2,1)=H(2)(3)(1,2)=0;
    H(3)(2)(2,2)=H(2)(3)(2,2)=0;
    H(3)(3)(2,1)=H(3)(3)(1,2)=0;
    H(3)(3)(2,2)=0;
    H(4)(1)(1,1)=H(1)(4)(1,1)=0;
    H(4)(1)(1,2)=H(1)(4)(2,1)=0;
    H(4)(2)(1,1)=H(2)(4)(1,1)=0;
    H(4)(2)(1,2)=H(2)(4)(2,1)=0;
    H(4)(3)(1,1)=H(3)(4)(1,1)=0;
    H(4)(3)(1,2)=H(3)(4)(2,1)=0;
    H(4)(4)(1,1)=0;
    H(4)(1)(2,1)=H(1)(4)(1,2)=0;
    H(4)(1)(2,2)=H(1)(4)(2,2)=0;
    H(4)(2)(2,1)=H(2)(4)(1,2)=0;
    H(4)(2)(2,2)=H(2)(4)(2,2)=0;
    H(4)(3)(2,1)=H(3)(4)(1,2)=0;
    H(4)(3)(2,2)=H(3)(4)(2,2)=0;
    H(4)(4)(2,1)=H(4)(4)(1,2)=0;
    H(4)(4)(2,2)=0;

//    LOG::cout<<"A "<<A<<std::endl;
    return A;
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1oo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
//    trap_cases.Append(2);
//    LOG::cout<<__FUNCTION__<<std::endl;
    T yvab=(a.y+b.y)/2,yba=b.y-a.y,yvcd=(c.y+d.y)/2,ydc=d.y-c.y,xba=b.x-a.x,xdc=d.x-c.x,xca=c.x-a.x;(void)yvab;(void)yba;(void)yvcd;(void)ydc;(void)xba;(void)xdc;(void)xca;

    T A=1./8*(xdc*xdc*xba*xba*ydc*ydc-4*xba*xba*yvcd*ydc*xdc*xdc-4*yvab*yba*xba*xba*xdc*xdc+4*xba*xba*xba*xba*ydc*ydc+yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xdc*xdc+4*ydc*ydc*xdc*xba*xba*xca-4*xdc*xba*xba*xba*ydc*yba-4*yba*xdc*ydc*xba*xca*xca+8*xba*xba*yba*xdc*ydc*xca+8*xba*yvab*yba*xdc*xdc*xca-8*yvcd*xdc*ydc*xba*xba*xca-8*yvab*xdc*xdc*yvcd*xba*xba-4*xdc*xdc*yba*yvcd*xba*xba+4*yvab*xdc*xdc*ydc*xba*xba-4*xba*yba*yba*xdc*xdc*xca+8*xba*xba*xba*xdc*yvcd*ydc+4*xba*xba*yvcd*yvcd*xdc*xdc-4*xdc*xba*xba*xba*ydc*ydc-8*xba*xba*xba*ydc*ydc*xca+4*xba*xba*ydc*ydc*xca*xca+2*yba*xba*xba*xdc*xdc*ydc+4*yba*yba*xdc*xdc*xca*xca)/xdc/xba/(-yba*xdc+ydc*xba);
    G(1)(1)=1./8*yba*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(4*ydc*xba*xba+2*yvcd*xba*xdc-ydc*xba*xdc-4*ydc*xba*xca-2*xba*yvab*xdc-3*xba*yba*xdc+2*yba*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    G(1)(2)=-1./8*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(4*ydc*xba*xba+2*yvcd*xba*xdc-ydc*xba*xdc-4*ydc*xba*xca-2*xba*yvab*xdc-3*xba*yba*xdc+2*yba*xdc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(2)(1)=1./8*(-4*xba*xba*xba*xba*ydc*ydc*ydc*xdc+8*xba*xba*xba*xba*xba*ydc*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc*xdc+4*yba*yba*xba*xba*xdc*xdc*xdc*yvcd-2*yba*yba*xba*xba*xdc*xdc*xdc*ydc+8*yba*ydc*ydc*xba*xba*xba*xdc*xdc-16*yba*yvcd*xba*xba*xba*xdc*xdc*ydc+8*yba*yba*xba*xba*xba*xdc*xdc*ydc-16*ydc*ydc*xba*xba*xba*xba*yba*xdc+4*yba*yvcd*xba*xba*ydc*xdc*xdc*xdc-4*yba*yvcd*yvcd*xba*xba*xdc*xdc*xdc-yba*ydc*ydc*xba*xba*xdc*xdc*xdc+8*yvab*xba*xba*yba*xdc*xdc*xdc*yvcd-4*yvab*xba*xba*yba*xdc*xdc*xdc*ydc-4*yvab*yvab*yba*xba*xba*xdc*xdc*xdc+4*yvab*yba*yba*xba*xba*xdc*xdc*xdc-yba*yba*yba*xba*xba*xdc*xdc*xdc-8*yvab*xba*xba*yba*xdc*xdc*ydc*xca+8*yba*yvcd*xba*xba*ydc*xdc*xdc*xca-4*yba*yba*xba*xba*xdc*xdc*ydc*xca+16*yba*ydc*ydc*xba*xba*xba*xdc*xca-4*yba*ydc*ydc*xba*xba*xdc*xdc*xca-8*xba*xba*xba*xba*ydc*ydc*ydc*xca+4*yba*yba*yba*xdc*xdc*xdc*xca*xca-8*yba*yba*xdc*xdc*ydc*xba*xca*xca)/xdc/xba/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(2)(2)=-1./8*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yvcd*xba*xdc+ydc*xba*xdc+4*ydc*xba*xca+2*xba*yvab*xdc-xba*yba*xdc-2*yba*xdc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(3)(1)=1./8*(3*ydc*ydc*ydc*xba*xba*xba*xdc*xdc-8*xba*xba*xba*xba*ydc*ydc*ydc*xdc+4*yba*yba*yba*xba*xdc*xdc*xdc*xdc+4*xba*xba*xba*xba*xba*ydc*ydc*ydc-4*yvcd*yvcd*ydc*xba*xba*xba*xdc*xdc-4*yvcd*ydc*ydc*xba*xba*xba*xdc*xdc-8*yvab*yba*yba*xba*xdc*xdc*xdc*xdc-12*yba*yba*xba*xba*xdc*xdc*xdc*ydc-4*yvab*ydc*ydc*xba*xba*xba*xdc*xdc+18*yba*ydc*ydc*xba*xba*xba*xdc*xdc+8*yvab*yvcd*xba*xba*xba*xdc*xdc*ydc-4*yba*yvcd*xba*xba*xba*xdc*xdc*ydc-4*yvab*yvab*xba*xba*xba*xdc*xdc*ydc+3*yba*yba*xba*xba*xba*xdc*xdc*ydc-8*ydc*ydc*xba*xba*xba*xba*yba*xdc+8*yba*yvcd*xba*xba*ydc*xdc*xdc*xdc-4*yba*ydc*ydc*xba*xba*xdc*xdc*xdc+8*yvab*xba*xba*yba*xdc*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc*xdc-8*yvab*xba*xba*yba*xdc*xdc*ydc*xca+8*yba*yvcd*xba*xba*ydc*xdc*xdc*xca+8*ydc*ydc*ydc*xba*xba*xba*xdc*xca-4*yba*yba*xba*xba*xdc*xdc*ydc*xca+16*yba*ydc*ydc*xba*xba*xba*xdc*xca-8*yba*ydc*ydc*xba*xba*xdc*xca*xca-20*yba*ydc*ydc*xba*xba*xdc*xdc*xca+16*yba*yba*xba*xdc*xdc*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xdc*xdc*xca+4*ydc*ydc*ydc*xba*xba*xba*xca*xca-8*xba*xba*xba*xba*ydc*ydc*ydc*xca)/xdc/xdc/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(3)(2)=-1./8*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(2*ydc*xba*xba-2*yvcd*xba*xdc-3*ydc*xba*xdc-2*ydc*xba*xca+2*xba*yvab*xdc-3*xba*yba*xdc+4*yba*xdc*xca+4*xdc*xdc*yba)/xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(4)(1)=-1./8*ydc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(2*ydc*xba*xba-2*yvcd*xba*xdc+ydc*xba*xdc-2*ydc*xba*xca+2*xba*yvab*xdc-3*xba*yba*xdc+4*yba*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xdc/xdc;
    G(4)(2)=1./8*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(2*ydc*xba*xba-2*yvcd*xba*xdc+ydc*xba*xdc-2*ydc*xba*xca+2*xba*yvab*xdc-3*xba*yba*xdc+4*yba*xdc*xca)/xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(1)(1)(1,1)=1./4*yba*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc+12*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-20*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca+4*ydc*xba*xba*xba*yvab*yba*xdc+24*yba*yba*xba*xba*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xba*xca+4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-4*xba*xba*xba*xba*ydc*ydc*ydc+4*yba*yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xba*ydc*xdc-11*yba*yba*xba*xba*xba*ydc*xdc-8*yvab*xba*xba*xba*xba*ydc*ydc+8*yba*xba*xba*xba*xba*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc+xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba/xba;
    H(1)(1)(2,1)=H(1)(1)(1,2)=-1./8*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc+4*yvcd*xba*xba*xba*yba*xdc*ydc+16*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-28*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-8*yvab*yba*xba*xba*xdc*xdc*yvcd+4*yvab*yba*xba*xba*xdc*xdc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca-4*yba*yba*xba*xba*xdc*xdc*yvcd+2*yba*yba*xba*xba*xdc*xdc*ydc-4*ydc*xba*xba*xba*yvab*yba*xdc+28*yba*yba*xba*xba*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xba*xca+4*xdc*xdc*yba*yvcd*yvcd*xba*xba+xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc-2*ydc*ydc*xba*xba*xba*yba*xdc-4*xba*xba*xba*xba*ydc*ydc*ydc+5*yba*yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xba*ydc*xdc-15*yba*yba*xba*xba*xba*ydc*xdc+4*yvab*yvab*yba*xba*xba*xdc*xdc+4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*xba*xba*ydc*ydc+12*yba*xba*xba*xba*xba*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc+8*yvab*xba*xba*yba*xdc*ydc*xca-8*yba*yvcd*xba*xba*xdc*ydc*xca+4*yba*ydc*ydc*xba*xba*xdc*xca+xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(1)(1)(2,2)=1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(1)(1,1)=H(1)(2)(1,1)=-1./4*yba*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc+12*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-12*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca+12*yba*yba*xba*xba*xdc*ydc*xca-4*yba*yba*yba*xdc*xdc*xba*xca+4*ydc*ydc*xba*xba*xba*yvab*xdc-2*xba*xba*xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*xba*xba*ydc*ydc+2*yba*xba*xba*xba*xba*ydc*ydc+4*xba*xba*xba*xba*yvcd*ydc*ydc+xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba/xba;
    H(2)(1)(1,2)=H(1)(2)(2,1)=-1./8*xdc*(-4*yba*yba*yba*xdc*xdc*xca*xca+8*yvcd*xba*xba*xba*yvab*xdc*ydc-12*yvcd*xba*xba*xba*yba*xdc*ydc-16*yba*ydc*ydc*xba*xba*xca*xca-8*yvab*ydc*ydc*xba*xba*xba*xca+12*yba*ydc*ydc*xba*xba*xba*xca+8*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+8*yvab*yba*xba*xba*xdc*xdc*yvcd-4*yvab*yba*xba*xba*xdc*xdc*ydc+12*yba*yba*ydc*xba*xdc*xca*xca+4*yba*yba*xba*xba*xdc*xdc*yvcd-2*yba*yba*xba*xba*xdc*xdc*ydc+12*ydc*xba*xba*xba*yvab*yba*xdc-4*yba*yba*xba*xba*xdc*ydc*xca-4*xdc*xdc*yba*yvcd*yvcd*xba*xba-xdc*xdc*yba*ydc*ydc*xba*xba-4*ydc*ydc*xba*xba*xba*yvab*xdc+6*ydc*ydc*xba*xba*xba*yba*xdc+3*yba*yba*yba*xba*xba*xdc*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-5*yba*yba*xba*xba*xba*ydc*xdc-4*yvab*yvab*yba*xba*xba*xdc*xdc-4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*yba*xdc*ydc*xca+8*yba*yvcd*xba*xba*xdc*ydc*xca-4*yba*ydc*ydc*xba*xba*xdc*xca-xba*xba*xba*xdc*ydc*ydc*ydc-4*ydc*ydc*ydc*xba*xba*xba*xca+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(2)(2)(1,1)=1./4*(-12*yba*yba*yba*xdc*xdc*xdc*ydc*xba*xca*xca+4*xba*xba*xba*yvcd*ydc*xdc*xdc*xdc*yba*yba+4*xba*xba*xba*ydc*ydc*xdc*xdc*xdc*yvab*yba-2*xba*xba*xba*ydc*ydc*xdc*xdc*xdc*yba*yba-12*xba*xba*xba*xba*xba*yba*xdc*ydc*ydc*ydc+4*xba*xba*xba*xba*xba*xba*ydc*ydc*ydc*ydc-4*xba*xba*xba*ydc*ydc*xdc*xdc*yba*yba*xca+12*yba*yba*xba*xba*xdc*xdc*ydc*ydc*xca*xca+4*ydc*ydc*ydc*xba*xba*xba*xdc*xdc*yba*xca+4*yvcd*yvcd*ydc*xba*xba*xba*xdc*xdc*xdc*yba-4*yvcd*ydc*ydc*xba*xba*xba*xdc*xdc*xdc*yba+12*xba*xba*xba*xba*xdc*xdc*yba*yba*ydc*ydc+ydc*ydc*ydc*xba*xba*xba*xdc*xdc*xdc*yba-8*yvcd*ydc*ydc*xba*xba*xba*xdc*xdc*yba*xca+8*yvab*yba*xba*xba*xba*xdc*xdc*ydc*ydc*xca-8*xba*xba*xba*yvcd*ydc*xdc*xdc*xdc*yvab*yba+4*yba*yba*yba*yba*xdc*xdc*xdc*xdc*xca*xca-3*xba*xba*xba*yba*yba*yba*xdc*xdc*xdc*ydc+4*xba*xba*xba*yvab*yvab*yba*xdc*xdc*xdc*ydc-4*xba*xba*xba*yvab*yba*yba*xdc*xdc*xdc*ydc)/xdc/xba/xba/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(1)(2,1)=H(1)(2)(1,2)=1./8*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc+16*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-20*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-8*yvab*yba*xba*xba*xdc*xdc*yvcd+4*yvab*yba*xba*xba*xdc*xdc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca+4*yba*yba*xba*xba*xdc*xdc*yvcd-2*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+20*yba*yba*xba*xba*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xba*xca+4*xdc*xdc*yba*yvcd*yvcd*xba*xba+xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-4*xba*xba*xba*xba*ydc*ydc*ydc+yba*yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xba*ydc*xdc-3*yba*yba*xba*xba*xba*ydc*xdc+4*yvab*yvab*yba*xba*xba*xdc*xdc-4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*xba*xba*ydc*ydc+4*yba*xba*xba*xba*xba*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc+8*yvab*xba*xba*yba*xdc*ydc*xca-8*yba*yvcd*xba*xba*xdc*ydc*xca+4*yba*ydc*ydc*xba*xba*xdc*xca+xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(2)(1)(2,2)=H(1)(2)(2,2)=1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(2)(2,1)=H(2)(2)(1,2)=1./8*xdc*(-4*yba*yba*yba*xdc*xdc*xca*xca+8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-16*yba*ydc*ydc*xba*xba*xca*xca-8*yvab*ydc*ydc*xba*xba*xba*xca+4*yba*ydc*ydc*xba*xba*xba*xca+8*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+8*yvab*yba*xba*xba*xdc*xdc*yvcd-4*yvab*yba*xba*xba*xdc*xdc*ydc+12*yba*yba*ydc*xba*xdc*xca*xca-4*yba*yba*xba*xba*xdc*xdc*yvcd+2*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+4*yba*yba*xba*xba*xdc*ydc*xca-4*xdc*xdc*yba*yvcd*yvcd*xba*xba-xdc*xdc*yba*ydc*ydc*xba*xba-4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-yba*yba*yba*xba*xba*xdc*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*yvab*yba*xba*xba*xdc*xdc+4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*yba*xdc*ydc*xca+8*yba*yvcd*xba*xba*xdc*ydc*xca-4*yba*ydc*ydc*xba*xba*xdc*xca-xba*xba*xba*xdc*ydc*ydc*ydc-4*ydc*ydc*ydc*xba*xba*xba*xca+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(2)(2)(2,2)=1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(1,1)=H(1)(3)(1,1)=-1./4*yba*(-8*yvcd*xba*xba*xba*yvab*xdc*ydc+4*yba*ydc*ydc*xba*xba*xca*xca+4*yvab*ydc*ydc*xba*xba*xba*xca-6*yba*ydc*ydc*xba*xba*xba*xca-4*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*yvab*yba*xba*xba*xdc*xdc*ydc-10*yba*yba*xba*xba*xdc*xdc*ydc+2*yba*yba*xba*xba*xdc*ydc*xca+2*xdc*xdc*yba*ydc*ydc*xba*xba+6*ydc*ydc*xba*xba*xba*yba*xdc-2*xba*xba*xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*xba*xba*ydc*ydc+2*yba*xba*xba*xba*xba*ydc*ydc+4*xba*xba*xba*xba*yvcd*ydc*ydc+4*xba*yba*yba*yba*xdc*xdc*xdc+4*yvab*xba*xba*yba*xdc*ydc*xca-4*yba*yvcd*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*ydc*xba*xca-6*yba*ydc*ydc*xba*xba*xdc*xca-xba*xba*xba*xdc*ydc*ydc*ydc+2*ydc*ydc*ydc*xba*xba*xba*xca-4*yba*yba*yba*xdc*xdc*xdc*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(3)(1)(1,2)=H(1)(3)(2,1)=1./4*(-8*yvcd*xba*xba*xba*yvab*xdc*ydc+4*yba*ydc*ydc*xba*xba*xca*xca+4*yvab*ydc*ydc*xba*xba*xba*xca-6*yba*ydc*ydc*xba*xba*xba*xca-4*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*yvab*yba*xba*xba*xdc*xdc*ydc-10*yba*yba*xba*xba*xdc*xdc*ydc+2*yba*yba*xba*xba*xdc*ydc*xca+2*xdc*xdc*yba*ydc*ydc*xba*xba+6*ydc*ydc*xba*xba*xba*yba*xdc-2*xba*xba*xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*xba*xba*ydc*ydc+2*yba*xba*xba*xba*xba*ydc*ydc+4*xba*xba*xba*xba*yvcd*ydc*ydc+4*xba*yba*yba*yba*xdc*xdc*xdc+4*yvab*xba*xba*yba*xdc*ydc*xca-4*yba*yvcd*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*ydc*xba*xca-6*yba*ydc*ydc*xba*xba*xdc*xca-xba*xba*xba*xdc*ydc*ydc*ydc+2*ydc*ydc*ydc*xba*xba*xba*xca-4*yba*yba*yba*xdc*xdc*xdc*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(1,1)=H(2)(3)(1,1)=1./4*(4*xba*xba*yvab*yba*yba*xdc*xdc*xdc*ydc*xca-4*xba*xba*xba*xba*xba*ydc*ydc*ydc*ydc*xca+4*xba*xba*xba*yvcd*ydc*xdc*xdc*xdc*yba*yba+4*ydc*xba*xba*xdc*xdc*xdc*xdc*yvab*yba*yba+12*yba*yba*yba*xdc*xdc*xdc*xdc*ydc*xba*xca-4*xba*xba*yvcd*ydc*xdc*xdc*xdc*yba*yba*xca-12*xba*xba*xba*ydc*ydc*xdc*xdc*xdc*yba*yba+2*ydc*xba*xba*xdc*xdc*xdc*xdc*yba*yba*yba+12*xba*xba*xba*xba*ydc*ydc*ydc*xdc*xdc*yba-12*xba*xba*xba*xba*xba*yba*xdc*ydc*ydc*ydc+4*xba*xba*xba*xba*xba*xba*ydc*ydc*ydc*ydc-4*xdc*ydc*ydc*ydc*ydc*xba*xba*xba*xba*xba-14*xba*xba*xba*ydc*ydc*xdc*xdc*yba*yba*xca+2*xba*xba*yba*yba*yba*xdc*xdc*xdc*ydc*xca-6*xba*xba*ydc*ydc*xdc*xdc*xdc*yba*yba*xca+4*yba*yba*xba*xba*xdc*xdc*ydc*ydc*xca*xca+2*ydc*ydc*ydc*xba*xba*xba*xdc*xdc*yba*xca+12*xba*xba*xba*xba*yba*xdc*ydc*ydc*ydc*xca+4*yvcd*yvcd*ydc*xba*xba*xba*xdc*xdc*xdc*yba-4*yvcd*ydc*xba*xba*xdc*xdc*xdc*xdc*yba*yba-4*yba*yba*yba*yba*xdc*xdc*xdc*xdc*xdc*xca+12*xba*xba*xba*xba*xdc*xdc*yba*yba*ydc*ydc-ydc*ydc*ydc*xba*xba*xba*xdc*xdc*xdc*yba+2*ydc*ydc*xba*xba*xdc*xdc*xdc*xdc*yba*yba-4*yvcd*ydc*ydc*xba*xba*xba*xdc*xdc*yba*xca+4*yvab*yba*xba*xba*xba*xdc*xdc*ydc*ydc*xca-8*xba*xba*xba*yvcd*ydc*xdc*xdc*xdc*yvab*yba-3*xba*xba*xba*yba*yba*yba*xdc*xdc*xdc*ydc+4*xba*xba*xba*yvab*yvab*yba*xdc*xdc*xdc*ydc-4*xba*xba*xba*yvab*yba*yba*xdc*xdc*xdc*ydc)/xdc/xdc/xba/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(1,2)=H(2)(3)(2,1)=1./4*(8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-4*yba*ydc*ydc*xba*xba*xca*xca-4*yvab*ydc*ydc*xba*xba*xba*xca+2*yba*ydc*ydc*xba*xba*xba*xca+4*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*yvab*yba*xba*xba*xdc*xdc*ydc+2*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+2*yba*yba*xba*xba*xdc*ydc*xca-2*xdc*xdc*yba*ydc*ydc*xba*xba-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*yba*xdc*ydc*xca+4*yba*yvcd*xba*xba*xdc*ydc*xca-12*yba*yba*xdc*xdc*ydc*xba*xca+6*yba*ydc*ydc*xba*xba*xdc*xca+xba*xba*xba*xdc*ydc*ydc*ydc-2*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(1,1)=1./4*(8*xba*xba*yvab*yba*yba*xdc*xdc*xdc*ydc*xca+4*ydc*ydc*ydc*ydc*xba*xba*xba*xba*xca*xca-8*xba*xba*xba*xba*xba*ydc*ydc*ydc*ydc*xca+4*xba*xba*xba*yvcd*ydc*xdc*xdc*xdc*yba*yba-4*xba*xba*xba*ydc*ydc*xdc*xdc*xdc*yvab*yba+8*ydc*xba*xba*xdc*xdc*xdc*xdc*yvab*yba*yba-8*xba*xba*yvcd*ydc*xdc*xdc*xdc*yba*yba*xca-22*xba*xba*xba*ydc*ydc*xdc*xdc*xdc*yba*yba+4*ydc*xba*xba*xdc*xdc*xdc*xdc*yba*yba*yba+24*xba*xba*xba*xba*ydc*ydc*ydc*xdc*xdc*yba-12*xba*xba*xba*xba*xba*yba*xdc*ydc*ydc*ydc+4*xba*xba*xba*xba*xba*xba*ydc*ydc*ydc*ydc+4*yba*yba*yba*yba*xdc*xdc*xdc*xdc*xdc*xdc+4*xdc*xdc*ydc*ydc*ydc*ydc*xba*xba*xba*xba-8*xdc*ydc*ydc*ydc*ydc*xba*xba*xba*xba*xba-24*xba*xba*xba*ydc*ydc*xdc*xdc*yba*yba*xca+4*xba*xba*yba*yba*yba*xdc*xdc*xdc*ydc*xca+20*xba*xba*ydc*ydc*xdc*xdc*xdc*yba*yba*xca+12*yba*yba*xba*xba*xdc*xdc*ydc*ydc*xca*xca-24*ydc*ydc*ydc*xba*xba*xba*xdc*xdc*yba*xca-12*ydc*ydc*ydc*xba*xba*xba*xdc*yba*xca*xca+24*xba*xba*xba*xba*yba*xdc*ydc*ydc*ydc*xca-12*yba*yba*yba*xdc*xdc*xdc*xdc*xdc*ydc*xba+8*xdc*ydc*ydc*ydc*ydc*xba*xba*xba*xba*xca+4*yvcd*yvcd*ydc*xba*xba*xba*xdc*xdc*xdc*yba+4*yvcd*ydc*ydc*xba*xba*xba*xdc*xdc*xdc*yba-8*yvcd*ydc*xba*xba*xdc*xdc*xdc*xdc*yba*yba+12*xba*xba*xba*xba*xdc*xdc*yba*yba*ydc*ydc-15*ydc*ydc*ydc*xba*xba*xba*xdc*xdc*xdc*yba+20*ydc*ydc*xba*xba*xdc*xdc*xdc*xdc*yba*yba-8*xba*xba*xba*yvcd*ydc*xdc*xdc*xdc*yvab*yba-3*xba*xba*xba*yba*yba*yba*xdc*xdc*xdc*ydc+4*xba*xba*xba*yvab*yvab*yba*xdc*xdc*xdc*ydc-4*xba*xba*xba*yvab*yba*yba*xdc*xdc*xdc*ydc)/xdc/xdc/xdc/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(2,1)=H(1)(3)(1,2)=1./4*yba*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(2,2)=H(1)(3)(2,2)=-1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(2,1)=H(2)(3)(1,2)=-1./4*(4*yba*yba*xdc*xdc*ydc*xca*xca-4*yvcd*xdc*xdc*xdc*yba*yba*xca+6*ydc*yba*yba*xdc*xdc*xdc*xca+4*yvab*yba*yba*xdc*xdc*xdc*xca-8*yvab*xdc*xdc*xdc*yvcd*xba*yba+12*yba*yba*xba*xba*xdc*xdc*ydc+4*yvcd*yvcd*yba*xba*xdc*xdc*xdc-ydc*ydc*yba*xba*xdc*xdc*xdc+4*yba*yba*xdc*xdc*xdc*yvcd*xba-12*yba*yba*xdc*xdc*xdc*ydc*xba+12*xdc*xdc*yba*ydc*ydc*xba*xba-12*ydc*ydc*xba*xba*xba*yba*xdc+4*xba*xba*xba*xba*ydc*ydc*ydc-3*xba*yba*yba*yba*xdc*xdc*xdc-14*yba*yba*xdc*xdc*ydc*xba*xca-2*ydc*ydc*yba*xba*xdc*xdc*xca+12*yba*ydc*ydc*xba*xba*xdc*xca+2*xdc*xdc*xdc*xdc*yba*yba*yba+4*yvab*yba*xdc*xdc*ydc*xba*xca-4*yvcd*yba*xba*ydc*xdc*xdc*xca-4*xba*xba*xba*xdc*ydc*ydc*ydc+4*xdc*xdc*xdc*xdc*yvab*yba*yba-4*yvcd*xdc*xdc*xdc*xdc*yba*yba+2*ydc*xdc*xdc*xdc*xdc*yba*yba-4*ydc*ydc*ydc*xba*xba*xba*xca+2*yba*yba*yba*xdc*xdc*xdc*xca+4*xba*yvab*yvab*yba*xdc*xdc*xdc-4*xba*yvab*yba*yba*xdc*xdc*xdc)/xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(2,2)=H(2)(3)(2,2)=-1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(2,1)=H(3)(3)(1,2)=-1./8*xba*(-12*yba*ydc*ydc*xba*xdc*xca*xca+16*yba*yba*xdc*xdc*ydc*xca*xca-8*yvcd*xdc*xdc*xdc*yba*yba*xca+28*ydc*yba*yba*xdc*xdc*xdc*xca+8*yvab*yba*yba*xdc*xdc*xdc*xca-8*yvab*xdc*xdc*xdc*yvcd*xba*yba+4*yvab*xdc*xdc*xdc*ydc*xba*yba-4*yvcd*yba*xba*ydc*xdc*xdc*xdc-4*yvab*yba*xba*xba*xdc*xdc*ydc+4*yvab*yvab*xba*xba*xdc*xdc*ydc+13*yba*yba*xba*xba*xdc*xdc*ydc+4*yvcd*yvcd*yba*xba*xdc*xdc*xdc-15*ydc*ydc*yba*xba*xdc*xdc*xdc+4*yba*yba*xdc*xdc*xdc*yvcd*xba-26*yba*yba*xdc*xdc*xdc*ydc*xba-4*yvab*xdc*xdc*ydc*ydc*xba*xba+26*xdc*xdc*yba*ydc*ydc*xba*xba-12*ydc*ydc*xba*xba*xba*yba*xdc+8*xdc*xba*xba*ydc*ydc*ydc*xca+4*xba*xba*xba*xba*ydc*ydc*ydc+5*xba*xba*ydc*ydc*ydc*xdc*xdc-3*xba*yba*yba*yba*xdc*xdc*xdc-28*yba*yba*xdc*xdc*ydc*xba*xca-28*ydc*ydc*yba*xba*xdc*xdc*xca+24*yba*ydc*ydc*xba*xba*xdc*xca+4*xdc*xdc*xdc*xdc*yba*yba*yba+8*yvab*yba*xdc*xdc*ydc*xba*xca-8*yvcd*yba*xba*ydc*xdc*xdc*xca-8*xba*xba*xba*xdc*ydc*ydc*ydc+8*xdc*xdc*xdc*xdc*yvab*yba*yba-8*yvcd*xdc*xdc*xdc*xdc*yba*yba+12*ydc*xdc*xdc*xdc*xdc*yba*yba+4*xba*xba*ydc*ydc*ydc*xca*xca-8*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca+4*xba*xba*yvcd*yvcd*ydc*xdc*xdc+4*xba*xba*yvcd*ydc*ydc*xdc*xdc+4*xba*yvab*yvab*yba*xdc*xdc*xdc-4*xba*yvab*yba*yba*xdc*xdc*xdc-8*yvab*xdc*xdc*yvcd*xba*xba*ydc+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xdc/xdc;
    H(3)(3)(2,2)=1./4*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(1,1)=H(1)(4)(1,1)=1./4*ydc*yba*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(1,2)=H(1)(4)(2,1)=-1./4*ydc*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(1,1)=H(2)(4)(1,1)=-1./4*ydc*(4*yba*yba*xdc*xdc*ydc*xca*xca-4*yvcd*xdc*xdc*xdc*yba*yba*xca+2*ydc*yba*yba*xdc*xdc*xdc*xca+4*yvab*yba*yba*xdc*xdc*xdc*xca-8*yvab*xdc*xdc*xdc*yvcd*xba*yba+4*yvab*xdc*xdc*xdc*ydc*xba*yba-4*yvcd*yba*xba*ydc*xdc*xdc*xdc+12*yba*yba*xba*xba*xdc*xdc*ydc+4*yvcd*yvcd*yba*xba*xdc*xdc*xdc+ydc*ydc*yba*xba*xdc*xdc*xdc+4*yba*yba*xdc*xdc*xdc*yvcd*xba-2*yba*yba*xdc*xdc*xdc*ydc*xba-12*ydc*ydc*xba*xba*xba*yba*xdc+4*xba*xba*xba*xba*ydc*ydc*ydc-3*xba*yba*yba*yba*xdc*xdc*xdc-14*yba*yba*xdc*xdc*ydc*xba*xca+2*ydc*ydc*yba*xba*xdc*xdc*xca+12*yba*ydc*ydc*xba*xba*xdc*xca+4*yvab*yba*xdc*xdc*ydc*xba*xca-4*yvcd*yba*xba*ydc*xdc*xdc*xca-4*ydc*ydc*ydc*xba*xba*xba*xca+2*yba*yba*yba*xdc*xdc*xdc*xca+4*xba*yvab*yvab*yba*xdc*xdc*xdc-4*xba*yvab*yba*yba*xdc*xdc*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xdc/xdc;
    H(4)(2)(1,2)=H(2)(4)(2,1)=-1./4*ydc*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(1,1)=H(3)(4)(1,1)=-1./4*ydc*xba*(-12*yba*ydc*ydc*xba*xdc*xca*xca+12*yba*yba*xdc*xdc*ydc*xca*xca-8*yvcd*xdc*xdc*xdc*yba*yba*xca+12*ydc*yba*yba*xdc*xdc*xdc*xca+8*yvab*yba*yba*xdc*xdc*xdc*xca-8*yvab*xdc*xdc*xdc*yvcd*xba*yba+12*yba*yba*xba*xba*xdc*xdc*ydc+4*yvcd*yvcd*yba*xba*xdc*xdc*xdc-ydc*ydc*yba*xba*xdc*xdc*xdc+4*yba*yba*xdc*xdc*xdc*yvcd*xba-12*yba*yba*xdc*xdc*xdc*ydc*xba+12*xdc*xdc*yba*ydc*ydc*xba*xba-12*ydc*ydc*xba*xba*xba*yba*xdc+4*xdc*xba*xba*ydc*ydc*ydc*xca+4*xba*xba*xba*xba*ydc*ydc*ydc-3*xba*yba*yba*yba*xdc*xdc*xdc-24*yba*yba*xdc*xdc*ydc*xba*xca-12*ydc*ydc*yba*xba*xdc*xdc*xca+24*yba*ydc*ydc*xba*xba*xdc*xca+2*xdc*xdc*xdc*xdc*yba*yba*yba-4*xba*xba*xba*xdc*ydc*ydc*ydc+4*xdc*xdc*xdc*xdc*yvab*yba*yba-4*yvcd*xdc*xdc*xdc*xdc*yba*yba+2*ydc*xdc*xdc*xdc*xdc*yba*yba+4*xba*xba*ydc*ydc*ydc*xca*xca-8*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca+4*xba*yvab*yvab*yba*xdc*xdc*xdc-4*xba*yvab*yba*yba*xdc*xdc*xdc)/(-yba*xdc+ydc*xba)
        /(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xdc/xdc/xdc;
    H(4)(3)(1,2)=H(3)(4)(2,1)=1./8*xba*(-12*yba*ydc*ydc*xba*xdc*xca*xca+16*yba*yba*xdc*xdc*ydc*xca*xca-8*yvcd*xdc*xdc*xdc*yba*yba*xca+12*ydc*yba*yba*xdc*xdc*xdc*xca+8*yvab*yba*yba*xdc*xdc*xdc*xca-8*yvab*xdc*xdc*xdc*yvcd*xba*yba+12*yvab*xdc*xdc*xdc*ydc*xba*yba-12*yvcd*yba*xba*ydc*xdc*xdc*xdc-4*yvab*yba*xba*xba*xdc*xdc*ydc+4*yvab*yvab*xba*xba*xdc*xdc*ydc+13*yba*yba*xba*xba*xdc*xdc*ydc+4*yvcd*yvcd*yba*xba*xdc*xdc*xdc+5*ydc*ydc*yba*xba*xdc*xdc*xdc+4*yba*yba*xdc*xdc*xdc*yvcd*xba-6*yba*yba*xdc*xdc*xdc*ydc*xba-4*yvab*xdc*xdc*ydc*ydc*xba*xba+2*xdc*xdc*yba*ydc*ydc*xba*xba-12*ydc*ydc*xba*xba*xba*yba*xdc+4*xba*xba*xba*xba*ydc*ydc*ydc-3*xba*xba*ydc*ydc*ydc*xdc*xdc-3*xba*yba*yba*yba*xdc*xdc*xdc-28*yba*yba*xdc*xdc*ydc*xba*xca-4*ydc*ydc*yba*xba*xdc*xdc*xca+24*yba*ydc*ydc*xba*xba*xdc*xca+8*yvab*yba*xdc*xdc*ydc*xba*xca-8*yvcd*yba*xba*ydc*xdc*xdc*xca+4*xba*xba*ydc*ydc*ydc*xca*xca-8*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca+4*xba*xba*yvcd*yvcd*ydc*xdc*xdc+4*xba*xba*yvcd*ydc*ydc*xdc*xdc+4*xba*yvab*yvab*yba*xdc*xdc*xdc-4*xba*yvab*yba*yba*xdc*xdc*xdc-8*yvab*xdc*xdc*yvcd*xba*xba*ydc+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xdc/xdc;
    H(4)(4)(1,1)=1./4*ydc*xba*(-12*yba*ydc*ydc*xba*xdc*xca*xca+12*yba*yba*xdc*xdc*ydc*xca*xca-8*yvcd*xdc*xdc*xdc*yba*yba*xca+4*ydc*yba*yba*xdc*xdc*xdc*xca+8*yvab*yba*yba*xdc*xdc*xdc*xca-8*yvab*xdc*xdc*xdc*yvcd*xba*yba+4*yvab*xdc*xdc*xdc*ydc*xba*yba-4*yvcd*yba*xba*ydc*xdc*xdc*xdc+12*yba*yba*xba*xba*xdc*xdc*ydc+4*yvcd*yvcd*yba*xba*xdc*xdc*xdc+ydc*ydc*yba*xba*xdc*xdc*xdc+4*yba*yba*xdc*xdc*xdc*yvcd*xba-2*yba*yba*xdc*xdc*xdc*ydc*xba-12*ydc*ydc*xba*xba*xba*yba*xdc+4*xba*xba*xba*xba*ydc*ydc*ydc-3*xba*yba*yba*yba*xdc*xdc*xdc-24*yba*yba*xdc*xdc*ydc*xba*xca+24*yba*ydc*ydc*xba*xba*xdc*xca+4*xba*xba*ydc*ydc*ydc*xca*xca-8*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca+4*xba*yvab*yvab*yba*xdc*xdc*xdc-4*xba*yvab*yba*yba*xdc*xdc*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xdc/xdc/xdc;
    H(4)(1)(2,1)=H(1)(4)(1,2)=-1./4*yba*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(2,2)=H(1)(4)(2,2)=1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(2,1)=H(2)(4)(1,2)=1./4*(4*yba*yba*xdc*xdc*ydc*xca*xca-4*yvcd*xdc*xdc*xdc*yba*yba*xca+2*ydc*yba*yba*xdc*xdc*xdc*xca+4*yvab*yba*yba*xdc*xdc*xdc*xca-8*yvab*xdc*xdc*xdc*yvcd*xba*yba+4*yvab*xdc*xdc*xdc*ydc*xba*yba-4*yvcd*yba*xba*ydc*xdc*xdc*xdc+12*yba*yba*xba*xba*xdc*xdc*ydc+4*yvcd*yvcd*yba*xba*xdc*xdc*xdc+ydc*ydc*yba*xba*xdc*xdc*xdc+4*yba*yba*xdc*xdc*xdc*yvcd*xba-2*yba*yba*xdc*xdc*xdc*ydc*xba-12*ydc*ydc*xba*xba*xba*yba*xdc+4*xba*xba*xba*xba*ydc*ydc*ydc-3*xba*yba*yba*yba*xdc*xdc*xdc-14*yba*yba*xdc*xdc*ydc*xba*xca+2*ydc*ydc*yba*xba*xdc*xdc*xca+12*yba*ydc*ydc*xba*xba*xdc*xca+4*yvab*yba*xdc*xdc*ydc*xba*xca-4*yvcd*yba*xba*ydc*xdc*xdc*xca-4*ydc*ydc*ydc*xba*xba*xba*xca+2*yba*yba*yba*xdc*xdc*xdc*xca+4*xba*yvab*yvab*yba*xdc*xdc*xdc-4*xba*yvab*yba*yba*xdc*xdc*xdc)/xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(2,2)=H(2)(4)(2,2)=1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(2,1)=H(3)(4)(1,2)=1./8*xba*(-12*yba*ydc*ydc*xba*xdc*xca*xca+16*yba*yba*xdc*xdc*ydc*xca*xca-8*yvcd*xdc*xdc*xdc*yba*yba*xca+20*ydc*yba*yba*xdc*xdc*xdc*xca+8*yvab*yba*yba*xdc*xdc*xdc*xca-8*yvab*xdc*xdc*xdc*yvcd*xba*yba-4*yvab*xdc*xdc*xdc*ydc*xba*yba+4*yvcd*yba*xba*ydc*xdc*xdc*xdc-4*yvab*yba*xba*xba*xdc*xdc*ydc+4*yvab*yvab*xba*xba*xdc*xdc*ydc+13*yba*yba*xba*xba*xdc*xdc*ydc+4*yvcd*yvcd*yba*xba*xdc*xdc*xdc-3*ydc*ydc*yba*xba*xdc*xdc*xdc+4*yba*yba*xdc*xdc*xdc*yvcd*xba-22*yba*yba*xdc*xdc*xdc*ydc*xba+4*yvab*xdc*xdc*ydc*ydc*xba*xba+22*xdc*xdc*yba*ydc*ydc*xba*xba-12*ydc*ydc*xba*xba*xba*yba*xdc+8*xdc*xba*xba*ydc*ydc*ydc*xca+4*xba*xba*xba*xba*ydc*ydc*ydc+xba*xba*ydc*ydc*ydc*xdc*xdc-3*xba*yba*yba*yba*xdc*xdc*xdc-28*yba*yba*xdc*xdc*ydc*xba*xca-20*ydc*ydc*yba*xba*xdc*xdc*xca+24*yba*ydc*ydc*xba*xba*xdc*xca+4*xdc*xdc*xdc*xdc*yba*yba*yba+8*yvab*yba*xdc*xdc*ydc*xba*xca-8*yvcd*yba*xba*ydc*xdc*xdc*xca-8*xba*xba*xba*xdc*ydc*ydc*ydc+8*xdc*xdc*xdc*xdc*yvab*yba*yba-8*yvcd*xdc*xdc*xdc*xdc*yba*yba+4*ydc*xdc*xdc*xdc*xdc*yba*yba+4*xba*xba*ydc*ydc*ydc*xca*xca-8*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca+4*xba*xba*yvcd*yvcd*ydc*xdc*xdc-4*xba*xba*yvcd*ydc*ydc*xdc*xdc+4*xba*yvab*yvab*yba*xdc*xdc*xdc-4*xba*yvab*yba*yba*xdc*xdc*xdc-8*yvab*xdc*xdc*yvcd*xba*xba*ydc+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xdc/xdc;
    H(4)(3)(2,2)=H(3)(4)(2,2)=-1./4*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(2,1)=H(4)(4)(1,2)=-1./8*xba*(-12*yba*ydc*ydc*xba*xdc*xca*xca+16*yba*yba*xdc*xdc*ydc*xca*xca-8*yvcd*xdc*xdc*xdc*yba*yba*xca+4*ydc*yba*yba*xdc*xdc*xdc*xca+8*yvab*yba*yba*xdc*xdc*xdc*xca-8*yvab*xdc*xdc*xdc*yvcd*xba*yba+4*yvab*xdc*xdc*xdc*ydc*xba*yba-4*yvcd*yba*xba*ydc*xdc*xdc*xdc-4*yvab*yba*xba*xba*xdc*xdc*ydc+4*yvab*yvab*xba*xba*xdc*xdc*ydc+13*yba*yba*xba*xba*xdc*xdc*ydc+4*yvcd*yvcd*yba*xba*xdc*xdc*xdc+ydc*ydc*yba*xba*xdc*xdc*xdc+4*yba*yba*xdc*xdc*xdc*yvcd*xba-2*yba*yba*xdc*xdc*xdc*ydc*xba+4*yvab*xdc*xdc*ydc*ydc*xba*xba-2*xdc*xdc*yba*ydc*ydc*xba*xba-12*ydc*ydc*xba*xba*xba*yba*xdc+4*xba*xba*xba*xba*ydc*ydc*ydc+xba*xba*ydc*ydc*ydc*xdc*xdc-3*xba*yba*yba*yba*xdc*xdc*xdc-28*yba*yba*xdc*xdc*ydc*xba*xca+4*ydc*ydc*yba*xba*xdc*xdc*xca+24*yba*ydc*ydc*xba*xba*xdc*xca+8*yvab*yba*xdc*xdc*ydc*xba*xca-8*yvcd*yba*xba*ydc*xdc*xdc*xca+4*xba*xba*ydc*ydc*ydc*xca*xca-8*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca+4*xba*xba*yvcd*yvcd*ydc*xdc*xdc-4*xba*xba*yvcd*ydc*ydc*xdc*xdc+4*xba*yvab*yvab*yba*xdc*xdc*xdc-4*xba*yvab*yba*yba*xdc*xdc*xdc-8*yvab*xdc*xdc*yvcd*xba*xba*ydc+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xdc/xdc;
    H(4)(4)(2,2)=1./4*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);

//    LOG::cout<<"A "<<A<<std::endl;
    return A;
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1uu(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
//    trap_cases.Append(3);
//    LOG::cout<<__FUNCTION__<<std::endl;
    T yvab=(a.y+b.y)/2,yba=b.y-a.y,yvcd=(c.y+d.y)/2,ydc=d.y-c.y,xba=b.x-a.x,xdc=d.x-c.x,xca=c.x-a.x;(void)yvab;(void)yba;(void)yvcd;(void)ydc;(void)xba;(void)xdc;(void)xca;

    T A=1./8*(-4*yvcd*yvcd*xba*xdc-ydc*ydc*xba*xdc+2*yba*xdc*ydc*xba-4*xba*yvab*yvab*xdc-4*ydc*yvab*xba*xdc-4*yvcd*yba*xba*xdc+8*yba*yvcd*xdc*xca-4*ydc*yba*xdc*xca+4*ydc*yba*xba*xca+8*yvcd*yvab*xba*xdc-4*xba*yvab*yba*xdc-xba*yba*yba*xdc-8*ydc*yvab*xba*xca+4*yvcd*xba*xdc*ydc+8*yvab*xba*xba*ydc-4*yba*xca*xca*ydc)/(-yba*xdc+ydc*xba);
    G(1)(1)=-1./8*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*yba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(1)(2)=1./8*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(2)(1)=1./8*(4*yba*yba*xdc*xdc*yvcd-2*yba*yba*xdc*xdc*ydc+8*yvab*ydc*ydc*xba*xba+4*yvab*yvab*yba*xdc*xdc+4*yvab*yba*yba*xdc*xdc+yba*yba*yba*xdc*xdc-8*yvab*yba*xdc*xdc*yvcd+4*yvab*yba*xdc*xdc*ydc+4*xdc*xdc*yba*yvcd*yvcd+xdc*xdc*yba*ydc*ydc-4*xdc*xdc*yba*yvcd*ydc-16*yvab*yba*xdc*ydc*xba+4*yba*ydc*ydc*xdc*xca-8*yba*yvcd*ydc*xdc*xca+8*yvab*yba*xdc*ydc*xca-4*yba*yba*xdc*ydc*xca+4*yba*ydc*ydc*xca*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(2)(2)=1./8*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(2*yvab*xdc-3*yba*xdc-2*yvcd*xdc+ydc*xdc+2*ydc*xca+2*ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(3)(1)=1./8*(8*yba*yvcd*xba*xdc*ydc-8*yvab*yvcd*xba*xba*ydc+4*yba*yvcd*xba*xba*ydc-4*yba*ydc*ydc*xba*xdc-4*yba*yba*ydc*xba*xdc-8*yba*yba*xdc*xdc*yvcd+4*yba*yba*xdc*xdc*ydc-4*yvab*ydc*ydc*xba*xba+2*yba*ydc*ydc*xba*xba+4*xba*xba*yvcd*yvcd*ydc-4*xba*xba*yvcd*ydc*ydc+xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*ydc+yba*yba*xba*xba*ydc+8*yvab*yba*xdc*ydc*xba-4*ydc*xba*xba*yvab*yba+4*yba*yba*ydc*xca*xca-8*yba*yvcd*xba*ydc*xca+8*yvab*yba*ydc*xba*xca-4*yba*ydc*ydc*xba*xca+8*yba*yba*xdc*ydc*xca-4*yba*yba*ydc*xba*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(3)(2)=-1./8*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-4*yba*xdc-2*yba*xca+2*yvcd*xba+3*ydc*xba-2*xba*yvab+xba*yba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(4)(1)=-1./8*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*ydc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(4)(2)=1./8*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(1)(1)(1,1)=-1./4*ydc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*yba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(1)(1)(2,1)=H(1)(1)(1,2)=1./8*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(ydc*xba+yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(1)(1)(2,2)=-1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(1)(1,1)=H(1)(2)(1,1)=-1./4*ydc*yba*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(1)(1,2)=H(1)(2)(2,1)=1./8*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(2*ydc*yba*xdc*xca+xdc*xdc*yba*ydc+ydc*ydc*xba*xdc+2*ydc*ydc*xba*xca+2*xba*xba*ydc*ydc-5*yba*xdc*ydc*xba+2*ydc*yvab*xba*xdc+2*yvab*xdc*xdc*yba+xdc*xdc*yba*yba-2*yvcd*xba*xdc*ydc-2*xdc*xdc*yba*yvcd)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(2)(1,1)=-1./4*ydc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*yba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(1)(2,1)=H(1)(2)(1,2)=-1./8*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*ydc*yba*xdc*xca-xdc*xdc*yba*ydc-ydc*ydc*xba*xdc-2*ydc*ydc*xba*xca+2*xba*xba*ydc*ydc-3*yba*xdc*ydc*xba-2*ydc*yvab*xba*xdc-2*yvab*xdc*xdc*yba+3*xdc*xdc*yba*yba+2*yvcd*xba*xdc*ydc+2*xdc*xdc*yba*yvcd)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(1)(2,2)=H(1)(2)(2,2)=-1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(2)(2,1)=H(2)(2)(1,2)=1./8*(4*yvcd*yvcd*xdc*xdc*xdc*yba+4*yvcd*xdc*xdc*xdc*yba*yba-2*ydc*yba*yba*xdc*xdc*xdc+ydc*ydc*yba*xdc*xdc*xdc-4*ydc*yvcd*xdc*xdc*xdc*yba+4*ydc*yvab*yba*xdc*xdc*xdc-8*yvcd*xdc*xdc*xdc*yvab*yba+8*yvab*yba*xdc*xdc*ydc*xca-8*yvcd*yba*ydc*xdc*xdc*xca+8*ydc*ydc*xba*xdc*yvab*xca-8*xdc*yvcd*ydc*ydc*xba*xca-2*ydc*ydc*yba*xba*xdc*xdc+4*yba*ydc*ydc*xdc*xca*xca-4*yba*yba*xdc*xdc*ydc*xca+4*ydc*ydc*yba*xdc*xdc*xca+4*xba*yvab*yvab*xdc*xdc*ydc+4*ydc*ydc*yvab*xba*xdc*xdc+xba*ydc*ydc*ydc*xdc*xdc+4*xba*ydc*ydc*ydc*xca*xca+4*xdc*xba*ydc*ydc*ydc*xca-8*yvcd*yvab*xba*ydc*xdc*xdc+4*yvcd*yba*xba*ydc*xdc*xdc+4*xba*yvcd*yvcd*ydc*xdc*xdc-4*xba*yvcd*ydc*ydc*xdc*xdc-4*yba*ydc*ydc*xba*xdc*xca-12*yba*ydc*ydc*xba*xba*xdc+13*yba*yba*xdc*xdc*ydc*xba-3*yba*yba*yba*xdc*xdc*xdc+4*ydc*ydc*ydc*xba*xba*xba+4*yvab*yvab*yba*xdc*xdc*xdc-4*yvab*yba*yba*xdc*xdc*xdc-4*yvab*yba*xdc*xdc*ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(2)(2,2)=-1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(1,1)=H(1)(3)(1,1)=1./4*yba*ydc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(1,2)=H(1)(3)(2,1)=-1./4*ydc*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(1,1)=H(2)(3)(1,1)=1./4*yba*ydc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(1,2)=H(2)(3)(2,1)=-1./4*ydc*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(1,1)=-1./4*yba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*ydc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(2,1)=H(1)(3)(1,2)=-1./4*yba*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(2,2)=H(1)(3)(2,2)=1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(2,1)=H(2)(3)(1,2)=-1./4*yba*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(2,2)=H(2)(3)(2,2)=1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(2,1)=H(3)(3)(1,2)=1./8*(4*xba*xba*yba*yba*xdc*yvcd+8*xba*yvab*yba*yba*xdc*xca-8*xba*yba*yba*xdc*yvcd*xca-4*yvab*yba*yba*xba*xba*xdc+4*yvab*yvab*yba*xba*xba*xdc+yba*yba*yba*xba*xba*xdc-8*xba*yba*yba*xdc*xdc*yvcd+4*yba*yba*yba*xdc*xca*xca-8*yvab*xba*xba*xdc*yvcd*yba+4*xba*xba*yba*xdc*yvcd*yvcd+9*yba*ydc*ydc*xba*xba*xdc-12*yba*yba*xdc*xdc*ydc*xba-2*yba*yba*xba*xba*xdc*ydc+4*yvcd*yvcd*ydc*xba*xba*xba+4*yvcd*ydc*ydc*xba*xba*xba-4*yba*yvcd*xba*xba*xdc*ydc+8*yba*yba*yba*xdc*xdc*xdc-3*ydc*ydc*ydc*xba*xba*xba-4*yvab*ydc*ydc*xba*xba*xba+2*yba*ydc*ydc*xba*xba*xba-8*yvab*yvcd*xba*xba*xba*ydc+4*yba*yvcd*xba*xba*xba*ydc+4*yvab*yvab*xba*xba*xba*ydc+yba*yba*xba*xba*xba*ydc-4*yba*yba*yba*xdc*xdc*xba+4*yvab*xba*xba*yba*xdc*ydc+8*yvab*yba*yba*xdc*xdc*xba-4*ydc*xba*xba*xba*yvab*yba-4*yba*yba*xba*xba*ydc*xca-4*yba*ydc*ydc*xba*xba*xca+4*yba*yba*ydc*xba*xca*xca+8*yvab*xba*xba*yba*ydc*xca-8*yba*yvcd*xba*xba*ydc*xca+4*yba*yba*ydc*xba*xdc*xca+8*yba*yba*yba*xdc*xdc*xca-4*xba*yba*yba*yba*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(2,2)=-1./4*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(1,1)=H(1)(4)(1,1)=-1./4*ydc*yba*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(1,2)=H(1)(4)(2,1)=1./4*ydc*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(1,1)=H(2)(4)(1,1)=-1./4*yba*ydc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(1,2)=H(2)(4)(2,1)=1./4*ydc*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(1,1)=H(3)(4)(1,1)=1./4*yba*ydc*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(1,2)=H(3)(4)(2,1)=-1./8*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(3*xba*xba*ydc*ydc-5*yba*xdc*ydc*xba+2*yvcd*yba*xba*xdc-2*ydc*yba*xba*xca-2*yba*yba*xdc*xca-2*xba*yvab*yba*xdc+xba*yba*yba*xdc+xba*xba*yba*ydc-2*yvab*xba*xba*ydc+2*xba*xba*yvcd*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(1,1)=-1./4*yba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*ydc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(2,1)=H(1)(4)(1,2)=1./4*yba*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(2,2)=H(1)(4)(2,2)=-1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(2,1)=H(2)(4)(1,2)=1./4*yba*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(2,2)=H(2)(4)(2,2)=-1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(2,1)=H(3)(4)(1,2)=-1./8*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-xba*xba*ydc*ydc+3*yba*xdc*ydc*xba+2*yvcd*yba*xba*xdc-2*ydc*yba*xba*xca-2*yba*yba*xdc*xca-4*xdc*xdc*yba*yba-2*xba*yvab*yba*xdc+xba*yba*yba*xdc+xba*xba*yba*ydc-2*yvab*xba*xba*ydc+2*xba*xba*yvcd*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(2,2)=H(3)(4)(2,2)=1./4*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(2,1)=H(4)(4)(1,2)=1./8*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(ydc*xba+yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(2,2)=-1./4*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);

//    LOG::cout<<"A "<<A<<std::endl;
    return A;
}

// Case 1 order: a c b d
template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1uo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
//    trap_cases.Append(4);
//    LOG::cout<<__FUNCTION__<<std::endl;
    T yvab=(a.y+b.y)/2,yba=b.y-a.y,yvcd=(c.y+d.y)/2,ydc=d.y-c.y,xba=b.x-a.x,xdc=d.x-c.x,xca=c.x-a.x;(void)yvab;(void)yba;(void)yvcd;(void)ydc;(void)xba;(void)xdc;(void)xca;
    T A=(T).5*(xba-xca)*(ydc*xba+2*yvcd*xdc-ydc*xdc-ydc*xca)/xdc;
    G(1)(1)=0;
    G(1)(2)=0;
    G(2)(1)=(T).5*(2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)/xdc;
    G(2)(2)=0;
    G(3)(1)=(T).5*(ydc*xba*xba-2*ydc*xba*xdc-2*ydc*xba*xca-2*yvcd*xdc*xdc+ydc*xdc*xdc+2*ydc*xdc*xca+ydc*xca*xca)/xdc/xdc;
    G(3)(2)=-(T).5*(xba-xca)*(xba-2*xdc-xca)/xdc;
    G(4)(1)=-(T).5*(xba-xca)*(xba-xca)*ydc/xdc/xdc;
    G(4)(2)=(T).5*(xba-xca)*(xba-xca)/xdc;
    H(1)(1)(1,1)=0;
    H(1)(1)(2,1)=H(1)(1)(1,2)=0;
    H(1)(1)(2,2)=0;
    H(2)(1)(1,1)=H(1)(2)(1,1)=0;
    H(2)(1)(1,2)=H(1)(2)(2,1)=0;
    H(2)(2)(1,1)=ydc/xdc;
    H(2)(1)(2,1)=H(1)(2)(1,2)=0;
    H(2)(1)(2,2)=H(1)(2)(2,2)=0;
    H(2)(2)(2,1)=H(2)(2)(1,2)=0;
    H(2)(2)(2,2)=0;
    H(3)(1)(1,1)=H(1)(3)(1,1)=0;
    H(3)(1)(1,2)=H(1)(3)(2,1)=0;
    H(3)(2)(1,1)=H(2)(3)(1,1)=ydc*(-xdc-xca+xba)/xdc/xdc;
    H(3)(2)(1,2)=H(2)(3)(2,1)=0;
    H(3)(3)(1,1)=(-xdc-xca+xba)*(-xdc-xca+xba)*ydc/xdc/xdc/xdc;
    H(3)(1)(2,1)=H(1)(3)(1,2)=0;
    H(3)(1)(2,2)=H(1)(3)(2,2)=0;
    H(3)(2)(2,1)=H(2)(3)(1,2)=-(-xdc-xca+xba)/xdc;
    H(3)(2)(2,2)=H(2)(3)(2,2)=0;
    H(3)(3)(2,1)=H(3)(3)(1,2)=-(T).5*(-2*xba*xdc-2*xba*xca+2*xdc*xdc+2*xdc*xca+xca*xca+xba*xba)/xdc/xdc;
    H(3)(3)(2,2)=0;
    H(4)(1)(1,1)=H(1)(4)(1,1)=0;
    H(4)(1)(1,2)=H(1)(4)(2,1)=0;
    H(4)(2)(1,1)=H(2)(4)(1,1)=-(xba-xca)*ydc/xdc/xdc;
    H(4)(2)(1,2)=H(2)(4)(2,1)=0;
    H(4)(3)(1,1)=H(3)(4)(1,1)=-ydc*(-xdc-xca+xba)*(xba-xca)/xdc/xdc/xdc;
    H(4)(3)(1,2)=H(3)(4)(2,1)=(T).5*(xba-xca)*(xba-xca)/xdc/xdc;
    H(4)(4)(1,1)=(xba-xca)*(xba-xca)*ydc/xdc/xdc/xdc;
    H(4)(1)(2,1)=H(1)(4)(1,2)=0;
    H(4)(1)(2,2)=H(1)(4)(2,2)=0;
    H(4)(2)(2,1)=H(2)(4)(1,2)=(xba-xca)/xdc;
    H(4)(2)(2,2)=H(2)(4)(2,2)=0;
    H(4)(3)(2,1)=H(3)(4)(1,2)=(T).5*(xba-xca)*(xba-2*xdc-xca)/xdc/xdc;
    H(4)(3)(2,2)=H(3)(4)(2,2)=0;
    H(4)(4)(2,1)=H(4)(4)(1,2)=-(T).5*(xba-xca)*(xba-xca)/xdc/xdc;
    H(4)(4)(2,2)=0;

//    LOG::cout<<"A "<<A<<std::endl;
    return A;
}

// Case 2 order: a c d b; ou = c over, d under
template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2ou(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
//    trap_cases.Append(5);
//    LOG::cout<<__FUNCTION__<<std::endl;
    T yvab=(a.y+b.y)/2,yba=b.y-a.y,yvcd=(c.y+d.y)/2,ydc=d.y-c.y,xba=b.x-a.x,xdc=d.x-c.x,xca=c.x-a.x;(void)yvab;(void)yba;(void)yvcd;(void)ydc;(void)xba;(void)xdc;(void)xca;

    T A=1./8*xdc*(-4*yba*yba*xba*xca+xba*xba*ydc*ydc-8*yvcd*yba*xba*xdc+yba*yba*xba*xba+4*ydc*yba*xba*xca+4*yvab*yvab*xba*xba+4*yba*yba*xca*xca+8*xba*yvab*yba*xca-2*xba*xba*yba*ydc-8*yvab*xba*xba*yvcd+4*yvab*xba*xba*ydc+4*xba*xba*yba*yvcd-4*yvab*yba*xba*xba+4*xba*xba*yvcd*yvcd+4*xba*xba*yvcd*ydc-8*yvcd*yba*xba*xca)/xba/(-yba*xdc+ydc*xba);
    G(1)(1)=1./8*yba*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(4*ydc*xba*xba+2*yvcd*xba*xdc-ydc*xba*xdc-4*ydc*xba*xca-2*xba*yvab*xdc-3*xba*yba*xdc+2*yba*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    G(1)(2)=-1./8*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(4*ydc*xba*xba+2*yvcd*xba*xdc-ydc*xba*xdc-4*ydc*xba*xca-2*xba*yvab*xdc-3*xba*yba*xdc+2*yba*xdc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(2)(1)=1./8*xdc*yba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yvcd*xba*xdc+ydc*xba*xdc+4*ydc*xba*xca+2*xba*yvab*xdc-xba*yba*xdc-2*yba*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    G(2)(2)=-1./8*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yvcd*xba*xdc+ydc*xba*xdc+4*ydc*xba*xca+2*xba*yvab*xdc-xba*yba*xdc-2*yba*xdc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(3)(1)=-1./8*(-4*yba*ydc*ydc*xba*xba*xdc+4*yba*yba*xdc*xdc*ydc*xba+4*yba*yba*xba*xba*xdc*ydc+4*yvcd*yvcd*ydc*xba*xba*xba+4*yvcd*ydc*ydc*xba*xba*xba-8*yba*yvcd*xba*xba*xdc*ydc+ydc*ydc*ydc*xba*xba*xba+4*yvab*ydc*ydc*xba*xba*xba-2*yba*ydc*ydc*xba*xba*xba-8*yvab*yvcd*xba*xba*xba*ydc+4*yba*yvcd*xba*xba*xba*ydc+4*yvab*yvab*xba*xba*xba*ydc+yba*yba*xba*xba*xba*ydc-4*yba*yba*yba*xdc*xdc*xba-8*yvab*xba*xba*yba*xdc*ydc+8*yvab*yba*yba*xdc*xdc*xba-4*ydc*xba*xba*xba*yvab*yba-4*yba*yba*xba*xba*ydc*xca+4*yba*ydc*ydc*xba*xba*xca+4*yba*yba*ydc*xba*xca*xca+8*yvab*xba*xba*yba*ydc*xca-8*yba*yvcd*xba*xba*ydc*xca-8*yba*yba*ydc*xba*xdc*xca+8*yba*yba*yba*xdc*xdc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(3)(2)=1./8*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(4)(1)=1./8*(-16*yba*yvcd*xba*xdc*ydc-8*yvab*yvcd*xba*xba*ydc+4*yba*yvcd*xba*xba*ydc+8*yba*yba*xdc*xdc*yvcd+4*yvab*ydc*ydc*xba*xba-2*yba*ydc*ydc*xba*xba+4*xba*xba*yvcd*yvcd*ydc+4*xba*xba*yvcd*ydc*ydc+xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*ydc+yba*yba*xba*xba*ydc-4*ydc*xba*xba*yvab*yba+4*yba*yba*ydc*xca*xca-8*yba*yvcd*xba*ydc*xca+8*yvab*yba*ydc*xba*xca+4*yba*ydc*ydc*xba*xca-4*yba*yba*ydc*xba*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(4)(2)=-1./8*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(2*yba*xdc-2*yba*xca+2*yvcd*xba-3*ydc*xba-2*xba*yvab+xba*yba)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(1)(1)(1,1)=1./4*yba*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc+12*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-20*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca+4*ydc*xba*xba*xba*yvab*yba*xdc+24*yba*yba*xba*xba*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xba*xca+4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-4*xba*xba*xba*xba*ydc*ydc*ydc+4*yba*yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xba*ydc*xdc-11*yba*yba*xba*xba*xba*ydc*xdc-8*yvab*xba*xba*xba*xba*ydc*ydc+8*yba*xba*xba*xba*xba*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc+xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba/xba;
    H(1)(1)(2,1)=H(1)(1)(1,2)=-1./8*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc+4*yvcd*xba*xba*xba*yba*xdc*ydc+16*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-28*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-8*yvab*yba*xba*xba*xdc*xdc*yvcd+4*yvab*yba*xba*xba*xdc*xdc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca-4*yba*yba*xba*xba*xdc*xdc*yvcd+2*yba*yba*xba*xba*xdc*xdc*ydc-4*ydc*xba*xba*xba*yvab*yba*xdc+28*yba*yba*xba*xba*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xba*xca+4*xdc*xdc*yba*yvcd*yvcd*xba*xba+xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc-2*ydc*ydc*xba*xba*xba*yba*xdc-4*xba*xba*xba*xba*ydc*ydc*ydc+5*yba*yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xba*ydc*xdc-15*yba*yba*xba*xba*xba*ydc*xdc+4*yvab*yvab*yba*xba*xba*xdc*xdc+4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*xba*xba*ydc*ydc+12*yba*xba*xba*xba*xba*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc+8*yvab*xba*xba*yba*xdc*ydc*xca-8*yba*yvcd*xba*xba*xdc*ydc*xca+4*yba*ydc*ydc*xba*xba*xdc*xca+xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(1)(1)(2,2)=1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(1)(1,1)=H(1)(2)(1,1)=-1./4*yba*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc+12*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-12*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca+12*yba*yba*xba*xba*xdc*ydc*xca-4*yba*yba*yba*xdc*xdc*xba*xca+4*ydc*ydc*xba*xba*xba*yvab*xdc-2*xba*xba*xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*xba*xba*ydc*ydc+2*yba*xba*xba*xba*xba*ydc*ydc+4*xba*xba*xba*xba*yvcd*ydc*ydc+xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba/xba;
    H(2)(1)(1,2)=H(1)(2)(2,1)=-1./8*xdc*(-4*yba*yba*yba*xdc*xdc*xca*xca+8*yvcd*xba*xba*xba*yvab*xdc*ydc-12*yvcd*xba*xba*xba*yba*xdc*ydc-16*yba*ydc*ydc*xba*xba*xca*xca-8*yvab*ydc*ydc*xba*xba*xba*xca+12*yba*ydc*ydc*xba*xba*xba*xca+8*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+8*yvab*yba*xba*xba*xdc*xdc*yvcd-4*yvab*yba*xba*xba*xdc*xdc*ydc+12*yba*yba*ydc*xba*xdc*xca*xca+4*yba*yba*xba*xba*xdc*xdc*yvcd-2*yba*yba*xba*xba*xdc*xdc*ydc+12*ydc*xba*xba*xba*yvab*yba*xdc-4*yba*yba*xba*xba*xdc*ydc*xca-4*xdc*xdc*yba*yvcd*yvcd*xba*xba-xdc*xdc*yba*ydc*ydc*xba*xba-4*ydc*ydc*xba*xba*xba*yvab*xdc+6*ydc*ydc*xba*xba*xba*yba*xdc+3*yba*yba*yba*xba*xba*xdc*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-5*yba*yba*xba*xba*xba*ydc*xdc-4*yvab*yvab*yba*xba*xba*xdc*xdc-4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*yba*xdc*ydc*xca+8*yba*yvcd*xba*xba*xdc*ydc*xca-4*yba*ydc*ydc*xba*xba*xdc*xca-xba*xba*xba*xdc*ydc*ydc*ydc-4*ydc*ydc*ydc*xba*xba*xba*xca+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(2)(2)(1,1)=-1./4*xdc*yba*(-4*yba*yba*yba*xdc*xdc*xca*xca+8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-12*yba*ydc*ydc*xba*xba*xca*xca-8*yvab*ydc*ydc*xba*xba*xba*xca+4*yba*ydc*ydc*xba*xba*xba*xca+8*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+12*yba*yba*ydc*xba*xdc*xca*xca+4*ydc*xba*xba*xba*yvab*yba*xdc-4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-xba*xba*xba*xdc*ydc*ydc*ydc-4*ydc*ydc*ydc*xba*xba*xba*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba/xba;
    H(2)(1)(2,1)=H(1)(2)(1,2)=1./8*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc+16*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-20*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-8*yvab*yba*xba*xba*xdc*xdc*yvcd+4*yvab*yba*xba*xba*xdc*xdc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca+4*yba*yba*xba*xba*xdc*xdc*yvcd-2*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+20*yba*yba*xba*xba*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xba*xca+4*xdc*xdc*yba*yvcd*yvcd*xba*xba+xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-4*xba*xba*xba*xba*ydc*ydc*ydc+yba*yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xba*ydc*xdc-3*yba*yba*xba*xba*xba*ydc*xdc+4*yvab*yvab*yba*xba*xba*xdc*xdc-4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*xba*xba*ydc*ydc+4*yba*xba*xba*xba*xba*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc+8*yvab*xba*xba*yba*xdc*ydc*xca-8*yba*yvcd*xba*xba*xdc*ydc*xca+4*yba*ydc*ydc*xba*xba*xdc*xca+xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(2)(1)(2,2)=H(1)(2)(2,2)=1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(2)(2,1)=H(2)(2)(1,2)=1./8*xdc*(-4*yba*yba*yba*xdc*xdc*xca*xca+8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-16*yba*ydc*ydc*xba*xba*xca*xca-8*yvab*ydc*ydc*xba*xba*xba*xca+4*yba*ydc*ydc*xba*xba*xba*xca+8*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+8*yvab*yba*xba*xba*xdc*xdc*yvcd-4*yvab*yba*xba*xba*xdc*xdc*ydc+12*yba*yba*ydc*xba*xdc*xca*xca-4*yba*yba*xba*xba*xdc*xdc*yvcd+2*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+4*yba*yba*xba*xba*xdc*ydc*xca-4*xdc*xdc*yba*yvcd*yvcd*xba*xba-xdc*xdc*yba*ydc*ydc*xba*xba-4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-yba*yba*yba*xba*xba*xdc*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*yvab*yba*xba*xba*xdc*xdc+4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*yba*xdc*ydc*xca+8*yba*yvcd*xba*xba*xdc*ydc*xca-4*yba*ydc*ydc*xba*xba*xdc*xca-xba*xba*xba*xdc*ydc*ydc*ydc-4*ydc*ydc*ydc*xba*xba*xba*xca+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(2)(2)(2,2)=1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(1,1)=H(1)(3)(1,1)=-1./4*yba*(-8*yvcd*xba*xba*xba*yvab*xdc*ydc+4*yba*ydc*ydc*xba*xba*xca*xca+4*yvab*ydc*ydc*xba*xba*xba*xca-6*yba*ydc*ydc*xba*xba*xba*xca-4*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*yvab*yba*xba*xba*xdc*xdc*ydc-10*yba*yba*xba*xba*xdc*xdc*ydc+2*yba*yba*xba*xba*xdc*ydc*xca+2*xdc*xdc*yba*ydc*ydc*xba*xba+6*ydc*ydc*xba*xba*xba*yba*xdc-2*xba*xba*xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*xba*xba*ydc*ydc+2*yba*xba*xba*xba*xba*ydc*ydc+4*xba*xba*xba*xba*yvcd*ydc*ydc+4*xba*yba*yba*yba*xdc*xdc*xdc+4*yvab*xba*xba*yba*xdc*ydc*xca-4*yba*yvcd*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*ydc*xba*xca-6*yba*ydc*ydc*xba*xba*xdc*xca-xba*xba*xba*xdc*ydc*ydc*ydc+2*ydc*ydc*ydc*xba*xba*xba*xca-4*yba*yba*yba*xdc*xdc*xdc*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(3)(1)(1,2)=H(1)(3)(2,1)=1./4*(-8*yvcd*xba*xba*xba*yvab*xdc*ydc+4*yba*ydc*ydc*xba*xba*xca*xca+4*yvab*ydc*ydc*xba*xba*xba*xca-6*yba*ydc*ydc*xba*xba*xba*xca-4*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*yvab*yba*xba*xba*xdc*xdc*ydc-10*yba*yba*xba*xba*xdc*xdc*ydc+2*yba*yba*xba*xba*xdc*ydc*xca+2*xdc*xdc*yba*ydc*ydc*xba*xba+6*ydc*ydc*xba*xba*xba*yba*xdc-2*xba*xba*xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*xba*xba*ydc*ydc+2*yba*xba*xba*xba*xba*ydc*ydc+4*xba*xba*xba*xba*yvcd*ydc*ydc+4*xba*yba*yba*yba*xdc*xdc*xdc+4*yvab*xba*xba*yba*xdc*ydc*xca-4*yba*yvcd*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*ydc*xba*xca-6*yba*ydc*ydc*xba*xba*xdc*xca-xba*xba*xba*xdc*ydc*ydc*ydc+2*ydc*ydc*ydc*xba*xba*xba*xca-4*yba*yba*yba*xdc*xdc*xdc*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(1,1)=H(2)(3)(1,1)=-1./4*yba*(8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-4*yba*ydc*ydc*xba*xba*xca*xca-4*yvab*ydc*ydc*xba*xba*xba*xca+2*yba*ydc*ydc*xba*xba*xba*xca+4*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*yvab*yba*xba*xba*xdc*xdc*ydc+2*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+2*yba*yba*xba*xba*xdc*ydc*xca-2*xdc*xdc*yba*ydc*ydc*xba*xba-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*yba*xdc*ydc*xca+4*yba*yvcd*xba*xba*xdc*ydc*xca-12*yba*yba*xdc*xdc*ydc*xba*xca+6*yba*ydc*ydc*xba*xba*xdc*xca+xba*xba*xba*xdc*ydc*ydc*ydc-2*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(3)(2)(1,2)=H(2)(3)(2,1)=1./4*(8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-4*yba*ydc*ydc*xba*xba*xca*xca-4*yvab*ydc*ydc*xba*xba*xba*xca+2*yba*ydc*ydc*xba*xba*xba*xca+4*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*yvab*yba*xba*xba*xdc*xdc*ydc+2*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+2*yba*yba*xba*xba*xdc*ydc*xca-2*xdc*xdc*yba*ydc*ydc*xba*xba-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*yba*xdc*ydc*xca+4*yba*yvcd*xba*xba*xdc*ydc*xca-12*yba*yba*xdc*xdc*ydc*xba*xca+6*yba*ydc*ydc*xba*xba*xdc*xca+xba*xba*xba*xdc*ydc*ydc*ydc-2*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(1,1)=1./4*yba*(8*yba*ydc*ydc*xba*xba*xdc-8*yba*yba*xdc*xdc*ydc*xba-4*yba*yba*xba*xba*xdc*ydc+4*yvcd*yvcd*ydc*xba*xba*xba+4*yvcd*ydc*ydc*xba*xba*xba-8*yba*yvcd*xba*xba*xdc*ydc+4*yba*yba*yba*xdc*xdc*xdc-3*ydc*ydc*ydc*xba*xba*xba-4*yvab*ydc*ydc*xba*xba*xba+2*yba*ydc*ydc*xba*xba*xba-8*yvab*yvcd*xba*xba*xba*ydc+4*yba*yvcd*xba*xba*xba*ydc+4*yvab*yvab*xba*xba*xba*ydc+yba*yba*xba*xba*xba*ydc+8*yvab*xba*xba*yba*xdc*ydc-4*ydc*xba*xba*xba*yvab*yba-4*yba*yba*xba*xba*ydc*xca-4*yba*ydc*ydc*xba*xba*xca+4*yba*yba*ydc*xba*xca*xca+8*yvab*xba*xba*yba*ydc*xca-8*yba*yvcd*xba*xba*ydc*xca+8*yba*yba*ydc*xba*xdc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(2,1)=H(1)(3)(1,2)=1./4*yba*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(2,2)=H(1)(3)(2,2)=-1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(2,1)=H(2)(3)(1,2)=1./4*yba*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(2,2)=H(2)(3)(2,2)=-1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(2,1)=H(3)(3)(1,2)=-1./8*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(ydc*xba+yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(2,2)=1./4*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(1,1)=H(1)(4)(1,1)=1./4*ydc*yba*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(1,2)=H(1)(4)(2,1)=-1./4*ydc*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(1,1)=H(2)(4)(1,1)=1./4*yba*ydc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(1,2)=H(2)(4)(2,1)=-1./4*ydc*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(1,1)=H(3)(4)(1,1)=-1./4*yba*ydc*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(1,2)=H(3)(4)(2,1)=1./8*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(xba*xba*ydc*ydc-5*yba*xdc*ydc*xba+2*yvcd*yba*xba*xdc-2*ydc*yba*xba*xca-2*yba*yba*xdc*xca+2*xdc*xdc*yba*yba-2*xba*yvab*yba*xdc+xba*yba*yba*xdc+xba*xba*yba*ydc-2*yvab*xba*xba*ydc+2*xba*xba*yvcd*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(1,1)=1./4*yba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*ydc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(2,1)=H(1)(4)(1,2)=-1./4*yba*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(2,2)=H(1)(4)(2,2)=1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(2,1)=H(2)(4)(1,2)=-1./4*yba*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(2,2)=H(2)(4)(2,2)=1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(2,1)=H(3)(4)(1,2)=1./8*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(-3*xba*xba*ydc*ydc+3*yba*xdc*ydc*xba+2*yvcd*yba*xba*xdc-2*ydc*yba*xba*xca-2*yba*yba*xdc*xca-2*xdc*xdc*yba*yba-2*xba*yvab*yba*xdc+xba*yba*yba*xdc+xba*xba*yba*ydc-2*yvab*xba*xba*ydc+2*xba*xba*yvcd*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(2,2)=H(3)(4)(2,2)=-1./4*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(2,1)=H(4)(4)(1,2)=-1./8*(4*xba*xba*yba*yba*xdc*yvcd+8*xba*yvab*yba*yba*xdc*xca-8*xba*yba*yba*xdc*yvcd*xca-4*yvab*yba*yba*xba*xba*xdc+4*yvab*yvab*yba*xba*xba*xdc+yba*yba*yba*xba*xba*xdc+4*yba*yba*yba*xdc*xca*xca-8*yvab*xba*xba*xdc*yvcd*yba+4*xba*xba*yba*xdc*yvcd*yvcd+13*yba*ydc*ydc*xba*xba*xdc-12*yba*yba*xdc*xdc*ydc*xba-2*yba*yba*xba*xba*xdc*ydc+4*yvcd*yvcd*ydc*xba*xba*xba-4*yvcd*ydc*ydc*xba*xba*xba-4*yba*yvcd*xba*xba*xdc*ydc+4*yba*yba*yba*xdc*xdc*xdc-3*ydc*ydc*ydc*xba*xba*xba+4*yvab*ydc*ydc*xba*xba*xba-2*yba*ydc*ydc*xba*xba*xba-8*yvab*yvcd*xba*xba*xba*ydc+4*yba*yvcd*xba*xba*xba*ydc+4*yvab*yvab*xba*xba*xba*ydc+yba*yba*xba*xba*xba*ydc+4*yvab*xba*xba*yba*xdc*ydc-4*ydc*xba*xba*xba*yvab*yba-4*yba*yba*xba*xba*ydc*xca+4*yba*ydc*ydc*xba*xba*xca+4*yba*yba*ydc*xba*xca*xca+8*yvab*xba*xba*yba*ydc*xca-8*yba*yvcd*xba*xba*ydc*xca+4*yba*yba*ydc*xba*xdc*xca-4*xba*yba*yba*yba*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(2,2)=1./4*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);

//    LOG::cout<<"A "<<A<<std::endl;
    return A;
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2oo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
//    trap_cases.Append(6);
////    LOG::cout<<__FUNCTION__<<std::endl;
    T yvab=(a.y+b.y)/2,yba=b.y-a.y,yvcd=(c.y+d.y)/2,ydc=d.y-c.y,xba=b.x-a.x,xdc=d.x-c.x,xca=c.x-a.x;(void)yvab;(void)yba;(void)yvcd;(void)ydc;(void)xba;(void)xdc;(void)xca;

    T A=-(T).5*xdc*(-2*yba*xca-2*xba*yvab+xba*yba-yba*xdc)/xba;
    G(1)(1)=-(T).5*xdc*yba*(-2*xca-xdc+2*xba)/xba/xba;
    G(1)(2)=(T).5*xdc*(-2*xca-xdc+2*xba)/xba;
    G(2)(1)=-(T).5*xdc*yba*(2*xca+xdc)/xba/xba;
    G(2)(2)=(T).5*xdc*(2*xca+xdc)/xba;
    G(3)(1)=(T).5*(-2*xba*yvab+xba*yba-2*yba*xca)/xba;
    G(3)(2)=0;
    G(4)(1)=-(T).5*(-2*xba*yvab+xba*yba-2*yba*xdc-2*yba*xca)/xba;
    G(4)(2)=0;
    H(1)(1)(1,1)=-xdc*yba*(-2*xca-xdc+2*xba)/xba/xba/xba;
    H(1)(1)(2,1)=H(1)(1)(1,2)=1./2*xdc*(-2*xca-xdc+2*xba)/xba/xba;
    H(1)(1)(2,2)=0;
    H(2)(1)(1,1)=H(1)(2)(1,1)=xdc*yba*(-xdc-2*xca+xba)/xba/xba/xba;
    H(2)(1)(1,2)=H(1)(2)(2,1)=1./2*xdc*(2*xca+xdc)/xba/xba;
    H(2)(2)(1,1)=xdc*yba*(2*xca+xdc)/xba/xba/xba;
    H(2)(1)(2,1)=H(1)(2)(1,2)=-1./2*xdc*(-2*xca-xdc+2*xba)/xba/xba;
    H(2)(1)(2,2)=H(1)(2)(2,2)=0;
    H(2)(2)(2,1)=H(2)(2)(1,2)=-1./2*xdc*(2*xca+xdc)/xba/xba;
    H(2)(2)(2,2)=0;
    H(3)(1)(1,1)=H(1)(3)(1,1)=(xba-xca)*yba/xba/xba;
    H(3)(1)(1,2)=H(1)(3)(2,1)=-(xba-xca)/xba;
    H(3)(2)(1,1)=H(2)(3)(1,1)=yba*xca/xba/xba;
    H(3)(2)(1,2)=H(2)(3)(2,1)=-xca/xba;
    H(3)(3)(1,1)=-yba/xba;
    H(3)(1)(2,1)=H(1)(3)(1,2)=0;
    H(3)(1)(2,2)=H(1)(3)(2,2)=0;
    H(3)(2)(2,1)=H(2)(3)(1,2)=0;
    H(3)(2)(2,2)=H(2)(3)(2,2)=0;
    H(3)(3)(2,1)=H(3)(3)(1,2)=0;
    H(3)(3)(2,2)=0;
    H(4)(1)(1,1)=H(1)(4)(1,1)=-yba*(-xdc-xca+xba)/xba/xba;
    H(4)(1)(1,2)=H(1)(4)(2,1)=(-xdc-xca+xba)/xba;
    H(4)(2)(1,1)=H(2)(4)(1,1)=-yba*(xdc+xca)/xba/xba;
    H(4)(2)(1,2)=H(2)(4)(2,1)=(xdc+xca)/xba;
    H(4)(3)(1,1)=H(3)(4)(1,1)=0;
    H(4)(3)(1,2)=H(3)(4)(2,1)=0;
    H(4)(4)(1,1)=yba/xba;
    H(4)(1)(2,1)=H(1)(4)(1,2)=0;
    H(4)(1)(2,2)=H(1)(4)(2,2)=0;
    H(4)(2)(2,1)=H(2)(4)(1,2)=0;
    H(4)(2)(2,2)=H(2)(4)(2,2)=0;
    H(4)(3)(2,1)=H(3)(4)(1,2)=0;
    H(4)(3)(2,2)=H(3)(4)(2,2)=0;
    H(4)(4)(2,1)=H(4)(4)(1,2)=0;
    H(4)(4)(2,2)=0;

//    LOG::cout<<"A "<<A<<std::endl;
    return A;
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2uu(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
//    trap_cases.Append(7);
//    LOG::cout<<__FUNCTION__<<std::endl;
    T yvab=(a.y+b.y)/2,yba=b.y-a.y,yvcd=(c.y+d.y)/2,ydc=d.y-c.y,xba=b.x-a.x,xdc=d.x-c.x,xca=c.x-a.x;(void)yvab;(void)yba;(void)yvcd;(void)ydc;(void)xba;(void)xdc;(void)xca;

    T A=yvcd*xdc;
    G(1)(1)=0;
    G(1)(2)=0;
    G(2)(1)=0;
    G(2)(2)=0;
    G(3)(1)=-yvcd;
    G(3)(2)=(T).5*xdc;
    G(4)(1)=yvcd;
    G(4)(2)=(T).5*xdc;
    H(1)(1)(1,1)=0;
    H(1)(1)(2,1)=H(1)(1)(1,2)=0;
    H(1)(1)(2,2)=0;
    H(2)(1)(1,1)=H(1)(2)(1,1)=0;
    H(2)(1)(1,2)=H(1)(2)(2,1)=0;
    H(2)(2)(1,1)=0;
    H(2)(1)(2,1)=H(1)(2)(1,2)=0;
    H(2)(1)(2,2)=H(1)(2)(2,2)=0;
    H(2)(2)(2,1)=H(2)(2)(1,2)=0;
    H(2)(2)(2,2)=0;
    H(3)(1)(1,1)=H(1)(3)(1,1)=0;
    H(3)(1)(1,2)=H(1)(3)(2,1)=0;
    H(3)(2)(1,1)=H(2)(3)(1,1)=0;
    H(3)(2)(1,2)=H(2)(3)(2,1)=0;
    H(3)(3)(1,1)=0;
    H(3)(1)(2,1)=H(1)(3)(1,2)=0;
    H(3)(1)(2,2)=H(1)(3)(2,2)=0;
    H(3)(2)(2,1)=H(2)(3)(1,2)=0;
    H(3)(2)(2,2)=H(2)(3)(2,2)=0;
    H(3)(3)(2,1)=H(3)(3)(1,2)=-(T).5;
    H(3)(3)(2,2)=0;
    H(4)(1)(1,1)=H(1)(4)(1,1)=0;
    H(4)(1)(1,2)=H(1)(4)(2,1)=0;
    H(4)(2)(1,1)=H(2)(4)(1,1)=0;
    H(4)(2)(1,2)=H(2)(4)(2,1)=0;
    H(4)(3)(1,1)=H(3)(4)(1,1)=0;
    H(4)(3)(1,2)=H(3)(4)(2,1)=(T).5;
    H(4)(4)(1,1)=0;
    H(4)(1)(2,1)=H(1)(4)(1,2)=0;
    H(4)(1)(2,2)=H(1)(4)(2,2)=0;
    H(4)(2)(2,1)=H(2)(4)(1,2)=0;
    H(4)(2)(2,2)=H(2)(4)(2,2)=0;
    H(4)(3)(2,1)=H(3)(4)(1,2)=-(T).5;
    H(4)(3)(2,2)=H(3)(4)(2,2)=0;
    H(4)(4)(2,1)=H(4)(4)(1,2)=(T).5;
    H(4)(4)(2,2)=0;

//    LOG::cout<<"A "<<A<<std::endl;
    return A;
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2uo(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
//    trap_cases.Append(8);
//    LOG::cout<<__FUNCTION__<<std::endl;
    T yvab=(a.y+b.y)/2,yba=b.y-a.y,yvcd=(c.y+d.y)/2,ydc=d.y-c.y,xba=b.x-a.x,xdc=d.x-c.x,xca=c.x-a.x;(void)yvab;(void)yba;(void)yvcd;(void)ydc;(void)xba;(void)xdc;(void)xca;

    T A=-1./8*xdc*(-4*yba*yba*xba*xca+xba*xba*ydc*ydc-4*yba*xdc*ydc*xba+8*yba*yba*xdc*xca+yba*yba*xba*xba-4*ydc*yba*xba*xca+4*yvab*yvab*xba*xba+4*yba*yba*xca*xca+8*xba*yvab*yba*xca+4*xdc*xdc*yba*yba+8*xba*yvab*yba*xdc-4*xba*yba*yba*xdc+2*xba*xba*yba*ydc-8*yvab*xba*xba*yvcd-4*yvab*xba*xba*ydc+4*xba*xba*yba*yvcd-4*yvab*yba*xba*xba+4*xba*xba*yvcd*yvcd-4*xba*xba*yvcd*ydc-8*yvcd*yba*xba*xca)/xba/(-yba*xdc+ydc*xba);
    G(1)(1)=-1./8*xdc*yba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(4*ydc*xba*xba+2*yvcd*xba*xdc-3*ydc*xba*xdc-4*ydc*xba*xca-2*xba*yvab*xdc-3*xba*yba*xdc+2*yba*xdc*xca+2*xdc*xdc*yba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    G(1)(2)=1./8*xdc*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(4*ydc*xba*xba+2*yvcd*xba*xdc-3*ydc*xba*xdc-4*ydc*xba*xca-2*xba*yvab*xdc-3*xba*yba*xdc+2*yba*xdc*xca+2*xdc*xdc*yba)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(2)(1)=-1./8*xdc*yba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(-2*yvcd*xba*xdc+3*ydc*xba*xdc+4*ydc*xba*xca+2*xba*yvab*xdc-xba*yba*xdc-2*yba*xdc*xca-2*xdc*xdc*yba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    G(2)(2)=1./8*xdc*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(-2*yvcd*xba*xdc+3*ydc*xba*xdc+4*ydc*xba*xca+2*xba*yvab*xdc-xba*yba*xdc-2*yba*xdc*xca-2*xdc*xdc*yba)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(3)(1)=1./8*(8*yba*yvcd*xba*xdc*ydc-8*yvab*yvcd*xba*xba*ydc+4*yba*yvcd*xba*xba*ydc-4*yba*ydc*ydc*xba*xdc-4*yba*yba*ydc*xba*xdc-8*yba*yba*xdc*xdc*yvcd+4*yba*yba*xdc*xdc*ydc-4*yvab*ydc*ydc*xba*xba+2*yba*ydc*ydc*xba*xba+4*xba*xba*yvcd*yvcd*ydc-4*xba*xba*yvcd*ydc*ydc+xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*ydc+yba*yba*xba*xba*ydc+8*yvab*yba*xdc*ydc*xba-4*ydc*xba*xba*yvab*yba+4*yba*yba*ydc*xca*xca-8*yba*yvcd*xba*ydc*xca+8*yvab*yba*ydc*xba*xca-4*yba*ydc*ydc*xba*xca+8*yba*yba*xdc*ydc*xca-4*yba*yba*ydc*xba*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(3)(2)=-1./8*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-4*yba*xdc-2*yba*xca+2*yvcd*xba+3*ydc*xba-2*xba*yvab+xba*yba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(4)(1)=-1./8*(-8*yba*ydc*ydc*xba*xba*xdc+16*yba*yba*xdc*xdc*ydc*xba-8*yba*yba*xba*xba*xdc*ydc+4*yvcd*yvcd*ydc*xba*xba*xba-4*yvcd*ydc*ydc*xba*xba*xba-8*yba*yba*yba*xdc*xdc*xdc+ydc*ydc*ydc*xba*xba*xba-4*yvab*ydc*ydc*xba*xba*xba+2*yba*ydc*ydc*xba*xba*xba-8*yvab*yvcd*xba*xba*xba*ydc+4*yba*yvcd*xba*xba*xba*ydc+4*yvab*yvab*xba*xba*xba*ydc+yba*yba*xba*xba*xba*ydc+4*yba*yba*yba*xdc*xdc*xba+16*yvab*xba*xba*yba*xdc*ydc-8*yvab*yba*yba*xdc*xdc*xba-4*ydc*xba*xba*xba*yvab*yba-4*yba*yba*xba*xba*ydc*xca-4*yba*ydc*ydc*xba*xba*xca+4*yba*yba*ydc*xba*xca*xca+8*yvab*xba*xba*yba*ydc*xca-8*yba*yvcd*xba*xba*ydc*xca+16*yba*yba*ydc*xba*xdc*xca-8*yba*yba*yba*xdc*xdc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    G(4)(2)=1./8*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(1)(1)(1,1)=-1./4*xdc*yba*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc+12*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-20*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca+24*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+24*yba*yba*xba*xba*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xba*xca-12*yba*yba*xdc*xdc*xdc*ydc*xba+12*xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc-22*ydc*ydc*xba*xba*xba*yba*xdc+4*xba*xba*xba*xba*ydc*ydc*ydc+4*yba*yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xba*ydc*xdc-11*yba*yba*xba*xba*xba*ydc*xdc-8*yvab*xba*xba*xba*xba*ydc*ydc+8*yba*xba*xba*xba*xba*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc-8*xba*yba*yba*yba*xdc*xdc*xdc-24*yba*yba*xdc*xdc*ydc*xba*xca+24*yba*ydc*ydc*xba*xba*xdc*xca+4*xdc*xdc*xdc*xdc*yba*yba*yba-3*xba*xba*xba*xdc*ydc*ydc*ydc-4*ydc*ydc*ydc*xba*xba*xba*xca+8*yba*yba*yba*xdc*xdc*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba/xba;
    H(1)(1)(2,1)=H(1)(1)(1,2)=1./8*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc+4*yvcd*xba*xba*xba*yba*xdc*ydc+16*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-28*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-8*yvab*yba*xba*xba*xdc*xdc*yvcd+4*yvab*yba*xba*xba*xdc*xdc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca-4*yba*yba*xba*xba*xdc*xdc*yvcd+26*yba*yba*xba*xba*xdc*xdc*ydc-4*ydc*xba*xba*xba*yvab*yba*xdc+28*yba*yba*xba*xba*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xba*xca-12*yba*yba*xdc*xdc*xdc*ydc*xba+4*xdc*xdc*yba*yvcd*yvcd*xba*xba+13*xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc-26*ydc*ydc*xba*xba*xba*yba*xdc+4*xba*xba*xba*xba*ydc*ydc*ydc+5*yba*yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xba*ydc*xdc-15*yba*yba*xba*xba*xba*ydc*xdc+4*yvab*yvab*yba*xba*xba*xdc*xdc+4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*xba*xba*ydc*ydc+12*yba*xba*xba*xba*xba*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc-8*xba*yba*yba*yba*xdc*xdc*xdc+8*yvab*xba*xba*yba*xdc*ydc*xca-8*yba*yvcd*xba*xba*xdc*ydc*xca-24*yba*yba*xdc*xdc*ydc*xba*xca+28*yba*ydc*ydc*xba*xba*xdc*xca+4*xdc*xdc*xdc*xdc*yba*yba*yba-3*xba*xba*xba*xdc*ydc*ydc*ydc-4*ydc*ydc*ydc*xba*xba*xba*xca+8*yba*yba*yba*xdc*xdc*xdc*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(1)(1)(2,2)=-1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(1)(1,1)=H(1)(2)(1,1)=1./4*xdc*yba*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc+12*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-12*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca+12*yba*yba*xba*xba*xdc*xdc*ydc+12*yba*yba*xba*xba*xdc*ydc*xca-4*yba*yba*yba*xdc*xdc*xba*xca-12*yba*yba*xdc*xdc*xdc*ydc*xba+12*xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc-12*ydc*ydc*xba*xba*xba*yba*xdc+2*xba*xba*xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*xba*xba*ydc*ydc+2*yba*xba*xba*xba*xba*ydc*ydc+4*xba*xba*xba*xba*yvcd*ydc*ydc-4*xba*yba*yba*yba*xdc*xdc*xdc-24*yba*yba*xdc*xdc*ydc*xba*xca+24*yba*ydc*ydc*xba*xba*xdc*xca+4*xdc*xdc*xdc*xdc*yba*yba*yba-3*xba*xba*xba*xdc*ydc*ydc*ydc-4*ydc*ydc*ydc*xba*xba*xba*xca+8*yba*yba*yba*xdc*xdc*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba/xba;
    H(2)(1)(1,2)=H(1)(2)(2,1)=1./8*xdc*(-4*yba*yba*yba*xdc*xdc*xca*xca+8*yvcd*xba*xba*xba*yvab*xdc*ydc-12*yvcd*xba*xba*xba*yba*xdc*ydc-16*yba*ydc*ydc*xba*xba*xca*xca-8*yvab*ydc*ydc*xba*xba*xba*xca+12*yba*ydc*ydc*xba*xba*xba*xca+8*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+8*yvab*yba*xba*xba*xdc*xdc*yvcd-4*yvab*yba*xba*xba*xdc*xdc*ydc+12*yba*yba*ydc*xba*xdc*xca*xca+4*yba*yba*xba*xba*xdc*xdc*yvcd-2*yba*yba*xba*xba*xdc*xdc*ydc+12*ydc*xba*xba*xba*yvab*yba*xdc-4*yba*yba*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*xdc*ydc*xba-4*xdc*xdc*yba*yvcd*yvcd*xba*xba-13*xdc*xdc*yba*ydc*ydc*xba*xba-4*ydc*ydc*xba*xba*xba*yvab*xdc+6*ydc*ydc*xba*xba*xba*yba*xdc+3*yba*yba*yba*xba*xba*xdc*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-5*yba*yba*xba*xba*xba*ydc*xdc-4*yvab*yvab*yba*xba*xba*xdc*xdc-4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*yba*xdc*ydc*xca+8*yba*yvcd*xba*xba*xdc*ydc*xca+24*yba*yba*xdc*xdc*ydc*xba*xca-28*yba*ydc*ydc*xba*xba*xdc*xca-4*xdc*xdc*xdc*xdc*yba*yba*yba+3*xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca-8*yba*yba*yba*xdc*xdc*xdc*xca+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(2)(2)(1,1)=1./4*xdc*yba*(-4*yba*yba*yba*xdc*xdc*xca*xca+8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-12*yba*ydc*ydc*xba*xba*xca*xca-8*yvab*ydc*ydc*xba*xba*xba*xca+4*yba*ydc*ydc*xba*xba*xba*xca+8*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+12*yba*yba*ydc*xba*xdc*xca*xca+4*ydc*xba*xba*xba*yvab*yba*xdc+12*yba*yba*xdc*xdc*xdc*ydc*xba-12*xdc*xdc*yba*ydc*ydc*xba*xba-4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc+24*yba*yba*xdc*xdc*ydc*xba*xca-24*yba*ydc*ydc*xba*xba*xdc*xca-4*xdc*xdc*xdc*xdc*yba*yba*yba+3*xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca-8*yba*yba*yba*xdc*xdc*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba/xba;
    H(2)(1)(2,1)=H(1)(2)(1,2)=-1./8*xdc*(4*yba*yba*yba*xdc*xdc*xca*xca-8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc+16*yba*ydc*ydc*xba*xba*xca*xca+8*yvab*ydc*ydc*xba*xba*xba*xca-20*yba*ydc*ydc*xba*xba*xba*xca-8*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc-8*yvab*yba*xba*xba*xdc*xdc*yvcd+4*yvab*yba*xba*xba*xdc*xdc*ydc-12*yba*yba*ydc*xba*xdc*xca*xca+4*yba*yba*xba*xba*xdc*xdc*yvcd+22*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+20*yba*yba*xba*xba*xdc*ydc*xca-8*yba*yba*yba*xdc*xdc*xba*xca-12*yba*yba*xdc*xdc*xdc*ydc*xba+4*xdc*xdc*yba*yvcd*yvcd*xba*xba+13*xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc-22*ydc*ydc*xba*xba*xba*yba*xdc+4*xba*xba*xba*xba*ydc*ydc*ydc+yba*yba*yba*xba*xba*xdc*xdc+4*yvab*yvab*xba*xba*xba*ydc*xdc-3*yba*yba*xba*xba*xba*ydc*xdc+4*yvab*yvab*yba*xba*xba*xdc*xdc-4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*xba*xba*ydc*ydc+4*yba*xba*xba*xba*xba*ydc*ydc+8*xba*xba*xba*xba*yvcd*ydc*ydc-8*xba*yba*yba*yba*xdc*xdc*xdc+8*yvab*xba*xba*yba*xdc*ydc*xca-8*yba*yvcd*xba*xba*xdc*ydc*xca-24*yba*yba*xdc*xdc*ydc*xba*xca+28*yba*ydc*ydc*xba*xba*xdc*xca+4*xdc*xdc*xdc*xdc*yba*yba*yba-3*xba*xba*xba*xdc*ydc*ydc*ydc-4*ydc*ydc*ydc*xba*xba*xba*xca+8*yba*yba*yba*xdc*xdc*xdc*xca-4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(2)(1)(2,2)=H(1)(2)(2,2)=-1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(2)(2)(2,1)=H(2)(2)(1,2)=-1./8*xdc*(-4*yba*yba*yba*xdc*xdc*xca*xca+8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-16*yba*ydc*ydc*xba*xba*xca*xca-8*yvab*ydc*ydc*xba*xba*xba*xca+4*yba*ydc*ydc*xba*xba*xba*xca+8*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+8*yvab*yba*xba*xba*xdc*xdc*yvcd-4*yvab*yba*xba*xba*xdc*xdc*ydc+12*yba*yba*ydc*xba*xdc*xca*xca-4*yba*yba*xba*xba*xdc*xdc*yvcd+2*yba*yba*xba*xba*xdc*xdc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+4*yba*yba*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*xdc*ydc*xba-4*xdc*xdc*yba*yvcd*yvcd*xba*xba-13*xdc*xdc*yba*ydc*ydc*xba*xba-4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-yba*yba*yba*xba*xba*xdc*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*yvab*yba*xba*xba*xdc*xdc+4*yvab*yba*yba*xba*xba*xdc*xdc-8*yvab*xba*xba*yba*xdc*ydc*xca+8*yba*yvcd*xba*xba*xdc*ydc*xca+24*yba*yba*xdc*xdc*ydc*xba*xca-28*yba*ydc*ydc*xba*xba*xdc*xca-4*xdc*xdc*xdc*xdc*yba*yba*yba+3*xba*xba*xba*xdc*ydc*ydc*ydc+4*ydc*ydc*ydc*xba*xba*xba*xca-8*yba*yba*yba*xdc*xdc*xdc*xca+4*xdc*xdc*yba*yvcd*xba*xba*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(2)(2)(2,2)=-1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(1,1)=H(1)(3)(1,1)=1./4*yba*ydc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(1,2)=H(1)(3)(2,1)=-1./4*ydc*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(1,1)=H(2)(3)(1,1)=1./4*yba*ydc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(1,2)=H(2)(3)(2,1)=-1./4*ydc*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(1,1)=-1./4*yba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*ydc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(2,1)=H(1)(3)(1,2)=-1./4*yba*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(1)(2,2)=H(1)(3)(2,2)=1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(2,1)=H(2)(3)(1,2)=-1./4*yba*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(2)(2,2)=H(2)(3)(2,2)=1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(2,1)=H(3)(3)(1,2)=1./8*(4*xba*xba*yba*yba*xdc*yvcd+8*xba*yvab*yba*yba*xdc*xca-8*xba*yba*yba*xdc*yvcd*xca-4*yvab*yba*yba*xba*xba*xdc+4*yvab*yvab*yba*xba*xba*xdc+yba*yba*yba*xba*xba*xdc-8*xba*yba*yba*xdc*xdc*yvcd+4*yba*yba*yba*xdc*xca*xca-8*yvab*xba*xba*xdc*yvcd*yba+4*xba*xba*yba*xdc*yvcd*yvcd+9*yba*ydc*ydc*xba*xba*xdc-12*yba*yba*xdc*xdc*ydc*xba-2*yba*yba*xba*xba*xdc*ydc+4*yvcd*yvcd*ydc*xba*xba*xba+4*yvcd*ydc*ydc*xba*xba*xba-4*yba*yvcd*xba*xba*xdc*ydc+8*yba*yba*yba*xdc*xdc*xdc-3*ydc*ydc*ydc*xba*xba*xba-4*yvab*ydc*ydc*xba*xba*xba+2*yba*ydc*ydc*xba*xba*xba-8*yvab*yvcd*xba*xba*xba*ydc+4*yba*yvcd*xba*xba*xba*ydc+4*yvab*yvab*xba*xba*xba*ydc+yba*yba*xba*xba*xba*ydc-4*yba*yba*yba*xdc*xdc*xba+4*yvab*xba*xba*yba*xdc*ydc+8*yvab*yba*yba*xdc*xdc*xba-4*ydc*xba*xba*xba*yvab*yba-4*yba*yba*xba*xba*ydc*xca-4*yba*ydc*ydc*xba*xba*xca+4*yba*yba*ydc*xba*xca*xca+8*yvab*xba*xba*yba*ydc*xca-8*yba*yvcd*xba*xba*ydc*xca+4*yba*yba*ydc*xba*xdc*xca+8*yba*yba*yba*xdc*xdc*xca-4*xba*yba*yba*yba*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(3)(3)(2,2)=-1./4*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(1,1)=H(1)(4)(1,1)=-1./4*yba*(-8*yvcd*xba*xba*xba*yvab*xdc*ydc+4*yba*ydc*ydc*xba*xba*xca*xca+4*yvab*ydc*ydc*xba*xba*xba*xca-6*yba*ydc*ydc*xba*xba*xba*xca-4*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc+12*yba*yba*xba*xba*xdc*xdc*ydc+2*yba*yba*xba*xba*xdc*ydc*xca-12*yba*yba*xdc*xdc*xdc*ydc*xba+12*xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc-12*ydc*ydc*xba*xba*xba*yba*xdc+2*xba*xba*xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*xba*xba*ydc*ydc+2*yba*xba*xba*xba*xba*ydc*ydc+4*xba*xba*xba*xba*yvcd*ydc*ydc-4*xba*yba*yba*yba*xdc*xdc*xdc+4*yvab*xba*xba*yba*xdc*ydc*xca-4*yba*yvcd*xba*xba*xdc*ydc*xca-12*yba*yba*xdc*xdc*ydc*xba*xca+14*yba*ydc*ydc*xba*xba*xdc*xca+4*xdc*xdc*xdc*xdc*yba*yba*yba-3*xba*xba*xba*xdc*ydc*ydc*ydc-2*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(4)(1)(1,2)=H(1)(4)(2,1)=1./4*(-8*yvcd*xba*xba*xba*yvab*xdc*ydc+4*yba*ydc*ydc*xba*xba*xca*xca+4*yvab*ydc*ydc*xba*xba*xba*xca-6*yba*ydc*ydc*xba*xba*xba*xca-4*yvcd*ydc*ydc*xba*xba*xba*xca+4*xba*xba*xba*xdc*yvcd*yvcd*ydc-4*xba*xba*xba*xdc*yvcd*ydc*ydc+12*yba*yba*xba*xba*xdc*xdc*ydc+2*yba*yba*xba*xba*xdc*ydc*xca-12*yba*yba*xdc*xdc*xdc*ydc*xba+12*xdc*xdc*yba*ydc*ydc*xba*xba+4*ydc*ydc*xba*xba*xba*yvab*xdc-12*ydc*ydc*xba*xba*xba*yba*xdc+2*xba*xba*xba*xba*ydc*ydc*ydc+4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*xba*xba*ydc*ydc+2*yba*xba*xba*xba*xba*ydc*ydc+4*xba*xba*xba*xba*yvcd*ydc*ydc-4*xba*yba*yba*yba*xdc*xdc*xdc+4*yvab*xba*xba*yba*xdc*ydc*xca-4*yba*yvcd*xba*xba*xdc*ydc*xca-12*yba*yba*xdc*xdc*ydc*xba*xca+14*yba*ydc*ydc*xba*xba*xdc*xca+4*xdc*xdc*xdc*xdc*yba*yba*yba-3*xba*xba*xba*xdc*ydc*ydc*ydc-2*ydc*ydc*ydc*xba*xba*xba*xca+4*yba*yba*yba*xdc*xdc*xdc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(1,1)=H(2)(4)(1,1)=-1./4*yba*(8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-4*yba*ydc*ydc*xba*xba*xca*xca-4*yvab*ydc*ydc*xba*xba*xba*xca+2*yba*ydc*ydc*xba*xba*xba*xca+4*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+2*yba*yba*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*xdc*ydc*xba-12*xdc*xdc*yba*ydc*ydc*xba*xba-4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*yba*xdc*ydc*xca+4*yba*yvcd*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*ydc*xba*xca-14*yba*ydc*ydc*xba*xba*xdc*xca-4*xdc*xdc*xdc*xdc*yba*yba*yba+3*xba*xba*xba*xdc*ydc*ydc*ydc+2*ydc*ydc*ydc*xba*xba*xba*xca-4*yba*yba*yba*xdc*xdc*xdc*xca)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/xba/xba;
    H(4)(2)(1,2)=H(2)(4)(2,1)=1./4*(8*yvcd*xba*xba*xba*yvab*xdc*ydc-4*yvcd*xba*xba*xba*yba*xdc*ydc-4*yba*ydc*ydc*xba*xba*xca*xca-4*yvab*ydc*ydc*xba*xba*xba*xca+2*yba*ydc*ydc*xba*xba*xba*xca+4*yvcd*ydc*ydc*xba*xba*xba*xca-4*xba*xba*xba*xdc*yvcd*yvcd*ydc+4*xba*xba*xba*xdc*yvcd*ydc*ydc+4*ydc*xba*xba*xba*yvab*yba*xdc+2*yba*yba*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*xdc*ydc*xba-12*xdc*xdc*yba*ydc*ydc*xba*xba-4*ydc*ydc*xba*xba*xba*yvab*xdc+2*ydc*ydc*xba*xba*xba*yba*xdc-4*yvab*yvab*xba*xba*xba*ydc*xdc-yba*yba*xba*xba*xba*ydc*xdc-4*yvab*xba*xba*yba*xdc*ydc*xca+4*yba*yvcd*xba*xba*xdc*ydc*xca+12*yba*yba*xdc*xdc*ydc*xba*xca-14*yba*ydc*ydc*xba*xba*xdc*xca-4*xdc*xdc*xdc*xdc*yba*yba*yba+3*xba*xba*xba*xdc*ydc*ydc*ydc+2*ydc*ydc*ydc*xba*xba*xba*xca-4*yba*yba*yba*xdc*xdc*xdc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(1,1)=H(3)(4)(1,1)=1./4*yba*ydc*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(1,2)=H(3)(4)(2,1)=-1./8*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(3*xba*xba*ydc*ydc-5*yba*xdc*ydc*xba+2*yvcd*yba*xba*xdc-2*ydc*yba*xba*xca-2*yba*yba*xdc*xca-2*xba*yvab*yba*xdc+xba*yba*yba*xdc+xba*xba*yba*ydc-2*yvab*xba*xba*ydc+2*xba*xba*yvcd*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(1,1)=-1./4*yba*(12*yba*ydc*ydc*xba*xba*xdc-12*yba*yba*xdc*xdc*ydc*xba+4*yvcd*yvcd*ydc*xba*xba*xba-4*yvcd*ydc*ydc*xba*xba*xba+4*yba*yba*yba*xdc*xdc*xdc-3*ydc*ydc*ydc*xba*xba*xba+4*yvab*ydc*ydc*xba*xba*xba-2*yba*ydc*ydc*xba*xba*xba-8*yvab*yvcd*xba*xba*xba*ydc+4*yba*yvcd*xba*xba*xba*ydc+4*yvab*yvab*xba*xba*xba*ydc+yba*yba*xba*xba*xba*ydc-4*ydc*xba*xba*xba*yvab*yba-4*yba*yba*xba*xba*ydc*xca+4*yba*ydc*ydc*xba*xba*xca+4*yba*yba*ydc*xba*xca*xca+8*yvab*xba*xba*yba*ydc*xca-8*yba*yvcd*xba*xba*ydc*xca)/xba/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(2,1)=H(1)(4)(1,2)=1./4*yba*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(1)(2,2)=H(1)(4)(2,2)=-1./4*xdc*(-2*yvab*xdc-yba*xdc+2*yvcd*xdc-ydc*xdc-2*ydc*xca+2*ydc*xba)*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(2,1)=H(2)(4)(1,2)=1./4*yba*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(2)(2,2)=H(2)(4)(2,2)=-1./4*xdc*(2*ydc*xca-2*yvcd*xdc+ydc*xdc+2*yvab*xdc-yba*xdc)*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(2,1)=H(3)(4)(1,2)=-1./8*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-xba*xba*ydc*ydc+3*yba*xdc*ydc*xba+2*yvcd*yba*xba*xdc-2*ydc*yba*xba*xca-2*yba*yba*xdc*xca-4*xdc*xdc*yba*yba-2*xba*yvab*yba*xdc+xba*yba*yba*xdc+xba*xba*yba*ydc-2*yvab*xba*xba*ydc+2*xba*xba*yvcd*ydc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(3)(2,2)=H(3)(4)(2,2)=1./4*xba*(ydc*xba-2*yba*xdc-2*yba*xca+2*yvcd*xba+xba*yba-2*xba*yvab)*xdc*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(2,1)=H(4)(4)(1,2)=1./8*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(ydc*xba+yba*xdc)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);
    H(4)(4)(2,2)=-1./4*xba*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*(-2*yba*xca-2*xba*yvab+xba*yba+2*yvcd*xba-ydc*xba)*xdc/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba)/(-yba*xdc+ydc*xba);

//    LOG::cout<<"A "<<A<<std::endl;
    return A;
}

// Case 1 order: a c b d; ou = c over, b under
template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_1(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    PHYSBAM_ASSERT(a.x<=c.x);
    PHYSBAM_ASSERT(c.x<=b.x);
    PHYSBAM_ASSERT(b.x<=d.x);
    bool co=(-a.y*b.x+a.x*b.y+c.y*b.x-c.y*a.x-b.y*c.x+a.y*c.x)>0;
    bool bo=(-d.y*b.x+c.y*b.x+b.y*d.x-b.y*c.x-c.y*d.x+c.x*d.y)>0;

    if(co){
        if(bo) return Trapezoid_Intersection_Area_Case_1oo(a,b,c,d,G,H);
        else return Trapezoid_Intersection_Area_Case_1ou(a,b,c,d,G,H);}
    else{
        if(bo) return Trapezoid_Intersection_Area_Case_1uo(a,b,c,d,G,H);
        else return Trapezoid_Intersection_Area_Case_1uu(a,b,c,d,G,H);}
}

// Case 2 order: a c d b; ou = c over, d under
template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area_Case_2(const TV& a,const TV& b,const TV& c,const TV& d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    PHYSBAM_ASSERT(a.x<=c.x);
    PHYSBAM_ASSERT(c.x<=d.x);
    PHYSBAM_ASSERT(d.x<=b.x);
    bool co=(-a.y*b.x+a.x*b.y+c.y*b.x-c.y*a.x-b.y*c.x+a.y*c.x)>0;
    bool Do=(-a.y*b.x+a.x*b.y+d.y*b.x-d.y*a.x-b.y*d.x+a.y*d.x)>0;

    if(co){
        if(Do) return Trapezoid_Intersection_Area_Case_2oo(a,b,c,d,G,H);
        else return Trapezoid_Intersection_Area_Case_2ou(a,b,c,d,G,H);}
    else{
        if(Do) return Trapezoid_Intersection_Area_Case_2uo(a,b,c,d,G,H);
        else return Trapezoid_Intersection_Area_Case_2uu(a,b,c,d,G,H);}
}

template<class T,class TV> T PhysBAM::
Trapezoid_Intersection_Area(TV a,TV b,TV c,TV d,VECTOR<TV,4>& G,VECTOR<VECTOR<MATRIX<T,2>,4>,4>& H)
{
    VECTOR<int,4> indices(1,2,3,4);
    T sign=1,A=0;
    if(a.x>b.x){sign=-sign;exchange(a,b);exchange(indices(1),indices(2));}
    if(c.x>d.x){sign=-sign;exchange(c,d);exchange(indices(3),indices(4));}
    if(a.x>c.x){exchange(a,c);exchange(b,d);exchange(indices(1),indices(3));exchange(indices(2),indices(4));}

    VECTOR<TV,4> tG;
    VECTOR<VECTOR<MATRIX<T,2>,4>,4> tH;

    if(b.x<c.x){G=tG;H=tH;return 0;}

    if(b.x<d.x) A=sign*Trapezoid_Intersection_Area_Case_1(a,b,c,d,tG,tH); // a c b d
    else A=sign*Trapezoid_Intersection_Area_Case_2(a,b,c,d,tG,tH); // a c d b
    for(int i=1;i<=4;i++){
        G(indices(i))=sign*tG(i);
        for(int j=1;j<=4;j++) H(indices(i))(indices(j))=sign*tH(i)(j);}

    return A;
}

template float Trapezoid_Intersection_Area<float,VECTOR<float,2> >(VECTOR<float,2>,VECTOR<float,2>,VECTOR<float,2>,VECTOR<float,2>,VECTOR<VECTOR<float,2>,4>&,VECTOR<VECTOR<MATRIX<float,2,2>,4>,4>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template double Trapezoid_Intersection_Area<double,VECTOR<double,2> >(VECTOR<double,2>,VECTOR<double,2>,VECTOR<double,2>,VECTOR<double,2>,VECTOR<VECTOR<double,2>,4>&,VECTOR<VECTOR<MATRIX<double,2,2>,4>,4>&);
#endif
