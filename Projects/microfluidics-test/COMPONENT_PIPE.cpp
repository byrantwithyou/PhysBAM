//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "COMPONENT_PIPE.h"
namespace PhysBAM{

//#####################################################################
// Function Make_Block
//#####################################################################
template<class T> CANONICAL_BLOCK<T>* COMPONENT_PIPE<T>::
Make_Block(const PIPE_KEY<T>& key)
{
    auto it=canonical_pipe_blocks.insert({key,{}});
    if(!it.second) return it.first->second;

    CANONICAL_BLOCK<T>* cb=new CANONICAL_BLOCK<T>;
    it.first->second=cb;
    
    int n=key.num_dofs;
    PHYSBAM_ASSERT(n%2);
    cb->X.Resize(2*n);
    cb->S.Resize(4*(n-1)+1);
    for(int i=0;i<n;i++)
    {
        T y=key.width/(n-1)*i-key.width/2;
        cb->X(i)=TV(0,y);
        cb->X(i+n)=TV(key.length,y);
    }

    for(int i=0;i<n-1;i++)
    {
        cb->E.Append({i,i+n,i+1});
        cb->E.Append({i+n,i+n+1,i+1});
        cb->S(i)={i,i+1};
        cb->S(i+(n-1))={i+n,i+n+1};
        cb->S(i+2*(n-1))={i+n,i+1};
        cb->S(i+3*(n-1))={i,i+n};
    }
    cb->S(4*(n-1))={n-1,2*n-1};

    cb->cross_sections.Append({{0,n},{0,n-1},false});
    cb->cross_sections.Append({{n,2*n},{n-1,2*(n-1)},true});
    cb->bc_v={0,n-1,n,2*n-1};
    cb->bc_e={3*(n-1),4*(n-1)};
    return it.first->second;
}
//#####################################################################
// Function Make_Component
//#####################################################################
template<class T> CANONICAL_COMPONENT<T>* COMPONENT_PIPE<T>::
Make_Component(const PIPE_KEY<T>& key)
{
    auto it=canonical_pipes.insert({key,{}});
    if(!it.second) return it.first->second;

    CANONICAL_COMPONENT<T>* cc=new CANONICAL_COMPONENT<T>;
    it.first->second=cc;

    T length=key.length;
    T offset=0;
    if(length>target_length*(T)1.5)
    {
        CANONICAL_BLOCK<T>* cb=Make_Block({key.num_dofs,key.width,target_length});
        while(length>target_length*(T)1.5)
        {
            cc->blocks.Append(
                {
                    cb,
                    {TV(offset,0)},
                    {{cc->blocks.m-1,CON_ID(1)},{cc->blocks.m+1,CON_ID(0)}}
                });
            length-=target_length;
            offset+=target_length;
        }
    }
    CANONICAL_BLOCK<T>* cb=Make_Block({key.num_dofs,key.width,length});
    cc->blocks.Append(
        {
            cb,
            {TV(offset,0)},
            {{cc->blocks.m-1,CON_ID(1)},{CC_BLOCK_ID(~1),CON_ID(0)}}
        });
    for(auto&bl:cc->blocks) bl.flags|=2;
    cc->blocks(CC_BLOCK_ID()).flags|=1;
    cc->blocks.Last().flags|=1;
    return cc;
}
template class COMPONENT_PIPE<float>;
template class COMPONENT_PIPE<double>;
}
