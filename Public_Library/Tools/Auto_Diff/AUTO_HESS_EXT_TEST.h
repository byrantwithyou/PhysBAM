//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AUTO_HESS_EXT_TEST
//##################################################################### 
#ifndef __AUTO_HESS_EXT_TEST__
#define __AUTO_HESS_EXT_TEST__

#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
namespace PhysBAM{
namespace HETERO_DIFF{

template<class TV,class VEC,class MAT,class EXTRA,int n>
bool Test(AUTO_HESS_EXT<TV,VEC,MAT>(*func)(const VECTOR<TV,n>&,EXTRA),EXTRA e)
{
    typedef typename TV::SCALAR T;
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
    T eps=1e-6;
    RANDOM_NUMBERS<T> rand;
    VECTOR<TV,n> x,dx;
    rand.Fill_Uniform(x,-1,1);
    rand.Fill_Uniform(dx,-eps,eps);
    VECTOR<TV,n> x2(x+dx);
    AUTO_HESS_EXT<TV,VEC,MAT> y=func(x,e);
    AUTO_HESS_EXT<TV,VEC,MAT> y2=func(x2,e);
    VECTOR<TV,n> g,g2;
    Extract(g,y.dx);
    Extract(g2,y2.dx);
    MATRIX<MATRIX<T,TV::m>,n> h,h2;
    Extract(h,y.ddx);
    Extract(h2,y2.ddx);

    T ta0=y2.x-y.x,tb0=0;
    for(int i=0;i<n;i++)
        tb0+=(g2(i)+g(i)).Dot(dx(i));
    tb0/=2;

    printf("first diff test %g (%g %g)\n",abs(ta0-tb0)/maxabs(ta0,tb0,1e-30),ta0,tb0);

    T ta1=0,tb1=0,tc1=0;
    for(int i=0;i<n;i++){
        TV r,s=g2(i)-g(i);
        for(int j=0;j<n;j++)
            r+=(h(i,j)*dx(j)+h2(i,j)*dx(j));
        r/=2;
        ta1+=r.Magnitude_Squared();
        tb1+=s.Magnitude_Squared();
        tc1+=(r-s).Magnitude_Squared();}
    ta1=sqrt(ta1);
    tb1=sqrt(tb1);
    tc1=sqrt(tc1);

    printf("second diff test %g (%g %g)\n",tc1/maxabs(ta1,tb1,1e-30),ta1,tb1);
    return true;
}
template<class TV,class VEC,class MAT,class EXTRA,int n>
bool Test(AUTO_HESS_EXT_VEC<TV,VEC,MAT>(*func)(const VECTOR<TV,n>&,EXTRA),EXTRA e)
{
    typedef typename TV::SCALAR T;
    RANDOM_NUMBERS<T> rand;
    TV v;
    rand.Fill_Uniform(v,-1,1);
    typedef TRIPLE<EXTRA,TV,AUTO_HESS_EXT_VEC<TV,VEC,MAT>(*)(const VECTOR<TV,n>&,EXTRA)> EX;
    EX pr(e,v,func);
    return Test(+[](const VECTOR<TV,n>& u,EX* p){return p->z(u,p->x).Dot(p->y);},&pr);
}
}
using HETERO_DIFF::Test;
}
#endif
