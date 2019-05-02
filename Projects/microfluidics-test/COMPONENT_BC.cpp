//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "COMPONENT_BC.h"
namespace PhysBAM{

//#####################################################################
// Function Make_Block
//#####################################################################
template<class T> auto COMPONENT_BC<T>::
Make_Block(int d,T w,T l,bool is_v) -> TRIP
{
    PIPE_KEY<T> key={d,w,l};
    auto it=canonical_bc_blocks[is_v].insert({key,{}});
    if(!it.second) return it.first->second;

    it.first->second.x=new CANONICAL_BLOCK<T>;
    auto* cb=it.first->second.x;

    int n=key.num_dofs;
    PHYSBAM_ASSERT(n%2);
    cb->X.Resize(2*n);
    cb->S.Resize(4*(n-1)+1);
    cb->ticks.Resize(cb->S.m,use_init,-7);
    for(int i=0;i<n;i++)
    {
        T y=key.width/(n-1)*i-key.width/2;
        cb->X(i)=TV(0,y);
        cb->X(i+n)=TV(key.length,y);
    }

    for(int i=0;i<n-1;i++)
    {
        if(i==0)
        {
            cb->E.Append({i,i+n,i+n+1});
            cb->E.Append({i+n+1,i+1,i});
        }
        else
        {
            cb->E.Append({i,i+n,i+1});
            cb->E.Append({i+n,i+n+1,i+1});
        }
        cb->S(i)={i,i+1};
        cb->S(i+(n-1))={i+n,i+n+1};
        if(i==0) cb->S(i+2*(n-1))={i+n+1,i};
        else cb->S(i+2*(n-1))={i+n,i+1};
        cb->S(i+3*(n-1))={i,i+n};
        cb->ticks(i)=0;
        cb->ticks(i+(n-1))=0;
        cb->ticks(i+2*(n-1))=0;
        cb->ticks(i+3*(n-1))=1;
    }
    cb->S(4*(n-1))={n-1,2*n-1};
    cb->ticks(4*(n-1))=1;

    cb->cross_sections.Append({{n,2*n},{n-1,2*(n-1)},true});
    cb->bc_v.Append(0);
    if(is_v) for(int i=1;i<n-1;i++) cb->bc_v.Append(i);
    cb->bc_v.Append(n-1);
    cb->bc_v.Append(n);
    cb->bc_v.Append(2*n-1);
    if(is_v) for(int i=0;i<n-1;i++) cb->bc_e.Append(i);
    cb->bc_e.Append(3*(n-1));
    cb->bc_e.Append(4*(n-1));
    it.first->second.y={0,n};
    it.first->second.z={0,n-1};
    cb->Compute_Element_Edges();
    return it.first->second;
}
template class COMPONENT_BC<double>;
}
