//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "COMPONENT_CHANGE.h"
namespace PhysBAM{

//#####################################################################
// Function Make_Canonical_Pipe_Change
//#####################################################################
template<class T> CANONICAL_COMPONENT<T>* COMPONENT_CHANGE<T>::
Make_Component(int d0,T w0,int d1,T w1,T l)
{
    PIPE_CHANGE_KEY<T> key={{d0,d1},{w0,w1},l};
    auto it=canonical_changes.insert({key,{}});
    if(!it.second) return it.first->second;

    int num_sec=rint(key.length/target_length);
    num_sec+=!num_sec;

    CANONICAL_COMPONENT<T>* cc=new CANONICAL_COMPONENT<T>;
    it.first->second=cc;

    T dw=w1-w0,wid=key.length/num_sec;
    int dd=d1-d0;
    for(int i=0;i<num_sec;i++)
    {
        T ow=w0+dw*i/num_sec;
        T nw=w0+dw*(i+1)/num_sec;
        int od=d0+dd*i/num_sec;
        int nd=d0+dd*(i+1)/num_sec;
        CANONICAL_BLOCK<T>* cb=Make_Block(od,ow,nd,nw,wid);
        cc->blocks.Append(
            {
                cb,
                {TV(i*wid,0)},
                {{cc->blocks.m-1,CON_ID(1)},{cc->blocks.m+1,CON_ID(0)}}
            });
    }
    cc->blocks.Last().connections(CON_ID(1)).id=CC_BLOCK_ID(~1);
    return cc;
}
template<class F> // func(a,b) means val(a)<val(b)
void Cross_Section_Topology(ARRAY<VECTOR<int,3> >& E,ARRAY<VECTOR<int,2> >& S,
    F func,int n0,int n1,ARRAY<int>& bc_e)
{
    for(int i=0;i<n0;i++) S.Append({i,i+1});
    for(int i=0;i<n1;i++) S.Append({i+n0,i+n0+1});
    int a=0,b=n0,c=n0+n1,ma=a+n0/2,mb=b+n1/2;
    bc_e.Append(S.m);
    while(a<n0 || b<c)
    {
        S.Append({a,b});
        // Ensure that the midpoints are topologically connected.
        bool need_a=(b>=c || (a==ma && b<mb));
        bool need_b=(a>=b || (b==mb && a<ma));
        if(need_a || (!need_b && func(a+1,b+1)))
        {
            E.Append({a,b,a+1});
            a++;
        }
        else
        {
            E.Append({a,b,b+1});
            b++;
        }
    }
    bc_e.Append(S.m);
    S.Append({a,b});
}
//#####################################################################
// Function Make_Canonical_Change_Block
//#####################################################################
template<class T> CANONICAL_BLOCK<T>* COMPONENT_CHANGE<T>::
Make_Block(int d0,T w0,int d1,T w1,T l)
{
    PIPE_CHANGE_KEY<T> key={{d0,d1},{w0,w1},l};
    auto it=canonical_change_blocks.insert({key,{}});
    if(!it.second) return it.first->second;

    it.first->second=new CANONICAL_BLOCK<T>;
    auto* cb=it.first->second;

    int n0=key.num_dofs[0];
    int n1=key.num_dofs[1];
    PHYSBAM_ASSERT(n0%2);
    PHYSBAM_ASSERT(n1%2);
    cb->X.Resize(n0+n1);
    for(int i=0;i<n0;i++)
        cb->X(i)=TV(0,key.width[0]/n0*i-key.width[0]/2);
    for(int i=0;i<n1;i++)
        cb->X(i+n0)=TV(key.length,key.width[1]/n1*i-key.width[1]/2);

    Cross_Section_Topology(cb->E,cb->S,
        [&cb](int a,int b){return cb->X(a+1).y<=cb->X(b+1).y;},n0,n1,cb->bc_e);
    cb->bc_v={0,n0-1,n0,n0+n1-1};

    cb->cross_sections.Append({{0,n0},{0,n0-1},false});
    cb->cross_sections.Append({{n0,2*n0},{n0-1,n0+n1-2},true});
    return it.first->second;
}
template struct COMPONENT_CHANGE<float>;
template struct COMPONENT_CHANGE<double>;
}
