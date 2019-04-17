//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VISITORS_FEM__
#define __VISITORS_FEM__

namespace PhysBAM{
// F(dof0,dof1,first_owns)
template<class F>
void Visit_Cross_Section_Dofs(INTERVAL<int> i0,INTERVAL<int> i1,bool uf0,bool uf1,bool m0,F func)
{
    int a=i1.min_corner,b=1,c=i0.min_corner,n=i0.Size();
    if(uf0==uf1)
    {
        a=i1.max_corner-1;
        b=-1;
    }
    if(uf0)
    {
        int m=n/2+(n&m0);
        for(int i=0;i<m;i++) func(c+i,a+b*i,true);
        for(int i=m;i<n;i++) func(c+i,a+b*i,false);
    }
    else
    {
        int m=n/2+(n&(1-m0));
        for(int i=0;i<m;i++) func(c+i,a+b*i,false);
        for(int i=m;i<n;i++) func(c+i,a+b*i,true);
    }
}

struct IRREGULAR_VISITOR
{
    BLOCK_ID b;
    int re,r0,r1;
    int ie,i0,i1;
    bool be,b0,b1; // regular owns e,v0,v1
    bool n0; // v0 not seen before for the irregular block
};

// F(const IRREGULAR_VISITOR& i)
template<class T,class F>
void Visit_Irregular_Cross_Section_Dofs(const ARRAY<BLOCK<T>,BLOCK_ID>& blocks,const IRREGULAR_CONNECTION& ic,F func)
{
    // edge: j=i*a+b
    // v0: j=i*a+c
    // v1: j=i*a+d
    auto& cs=blocks(ic.regular).block->cross_sections(ic.con_id);
    int a=1,b=0,c=0,d=1,n=ic.edge_on.m;
    if(!cs.own_first)
    {
        a=-1;
        b=n-1;
        d=n-1;
        c=n;
    }
    b+=cs.e.min_corner;
    c+=cs.v.min_corner;
    d+=cs.v.min_corner;

    int irreg_v0=n/2+1;
    int irreg_v1=n/2;
    BLOCK_ID last_b(-1);
    for(int i=0;i<n;i++)
    {
        bool b1=i<irreg_v1;
        bool b0=i<irreg_v0;
        auto& eo=ic.edge_on(i);
        IRREGULAR_VISITOR iv={eo.b,a*i+b,a*i+c,a*i+d,eo.e,eo.v0,eo.v1,b1,b0,b1,eo.b!=last_b};
        func(iv);
        last_b=eo.b;
    }
}
template<class CS,class FV,class FE>
void Visit_Regular_Cross_Section_Dofs(const CS& cs0,const CS& cs1,bool m0,FV func_v,FE func_e)
{
    Visit_Cross_Section_Dofs(cs0.v,cs1.v,cs0.own_first,cs1.own_first,m0,func_v);
    Visit_Cross_Section_Dofs(cs0.e,cs1.e,cs0.own_first,cs1.own_first,m0,func_e);
}
}
#endif
