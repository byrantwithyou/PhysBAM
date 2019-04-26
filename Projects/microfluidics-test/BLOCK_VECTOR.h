//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BLOCK_VECTOR__
#define __BLOCK_VECTOR__
#include <Core/Arrays/ARRAY.h>
#include "COMMON.h"

namespace PhysBAM{

template<class TV>
struct BLOCK_VECTOR
{
    typedef typename TV::SCALAR T;
    DOF_COUNTS n;
    ARRAY<T> V;

    void Resize()
    {
        V.Resize(TV::m*n.v+TV::m*n.e+n.p);
    }

    TV Get_v(int i) const {TV r;for(int j=0;j<TV::m;j++) r(j)=V(TV::m*i+j);return r;}
    TV Get_e(int i) const {return Get_v(i+n.v);}
    TV Get_u(int i,int e) const {return Get_v(i+e*n.v);}
    T Get_p(int i) const {return V(TV::m*(n.v+n.e)+i);}
    void Add_v(int i,TV u) {for(int j=0;j<TV::m;j++) V(TV::m*i+j)+=u(j);}
    void Add_e(int i,TV u) {Add_v(i+n.v,u);}
    void Add_u(int i,int e,TV u) {Add_v(i+e*n.v,u);}
    void Add_p(int i,T u) {V(TV::m*(n.v+n.e)+i)+=u;}
    void Set_v(int i,TV u) {for(int j=0;j<TV::m;j++) V(TV::m*i+j)=u(j);}
    void Set_e(int i,TV u) {Set_v(i+n.v,u);}
    void Set_u(int i,int e,TV u) {Set_v(i+e*n.v,u);}
    void Set_p(int i,T u) {V(TV::m*(n.v+n.e)+i)=u;}

    void Transform(const MATRIX<T,2>& M,T scale_p)
    {
        for(int i=0;i<n.v;i++) Set_v(i,M*Get_v(i));
        for(int i=0;i<n.e;i++) Set_e(i,M*Get_e(i));
        for(int i=0;i<n.p;i++) Set_p(i,scale_p*Get_p(i));
    }
};

}
#endif
