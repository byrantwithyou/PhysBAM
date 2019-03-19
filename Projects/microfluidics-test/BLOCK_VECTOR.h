//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BLOCK_VECTOR__
#define __BLOCK_VECTOR__
#include <Core/Arrays/ARRAY.h>
#include "COMMON.h"

namespace PhysBAM{

template<class T>
struct BLOCK_VECTOR
{
    typedef VECTOR<T,2> TV;
    DOF_COUNTS n;
    ARRAY<T> V;

    void Resize()
    {
        V.Resize(2*n.v+2*n.e+n.p);
    }

    TV Get_v(int i) const {return TV(V(2*i),V(2*i+1));}
    TV Get_e(int i) const {return Get_v(i+n.v);}
    TV Get_u(int i,int e) const {return Get_v(i+e*n.v);}
    T Get_p(int i) const {return V(2*(n.v+n.e)+i);}
    void Add_v(int i,TV u) {V(2*i)+=u(0);V(2*i+1)+=u(1);}
    void Add_e(int i,TV u) {Add_v(i+n.v,u);}
    void Add_u(int i,int e,TV u) {Add_v(i+e*n.v,u);}
    void Add_p(int i,T u) {V(2*(n.v+n.e)+i)+=u;}
    void Set_v(int i,TV u) {V(2*i)=u(0);V(2*i+1)=u(1);}
    void Set_e(int i,TV u) {Set_v(i+n.v,u);}
    void Set_u(int i,int e,TV u) {Set_v(i+e*n.v,u);}
    void Set_p(int i,T u) {V(2*(n.v+n.e)+i)=u;}
};

}
#endif
