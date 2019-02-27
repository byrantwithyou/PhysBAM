//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BLOCK_VECTOR__
#define __BLOCK_VECTOR__
#include <Core/Arrays/ARRAY.h>

namespace PhysBAM{

template<class T>
struct BLOCK_VECTOR
{
    typedef VECTOR<T,2> TV;
    int nv,ne,np;
    ARRAY<T> V;

    void Resize()
    {
        V.Resize(2*nv+2*ne+np);
    }

    TV Get_v(int i) const {return TV(V(2*i),V(2*i+1));}
    TV Get_e(int i) const {return Get_v(i+nv);}
    TV Get_u(int i,int e) const {return Get_v(i+e*nv);}
    T Get_p(int i) const {return V(2*(nv+ne)+i);}
    void Add_v(int i,TV u) {V(2*i)+=u(0);V(2*i+1)+=u(1);}
    void Add_e(int i,TV u) {Add_v(i+nv,u);}
    void Add_u(int i,int e,TV u) {Add_v(i+e*nv,u);}
    void Add_p(int i,T u) {V(2*(nv+ne)+i)+=u;}
};

}
#endif
