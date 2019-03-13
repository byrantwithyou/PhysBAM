//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPONENT_CHANGE__
#define __COMPONENT_CHANGE__

#include <map>
#include "CANONICAL_BLOCK.h"
#include "CANONICAL_COMPONENT.h"

namespace PhysBAM{

extern double comp_tol;

template<class T>
struct PIPE_CHANGE_KEY
{
    int num_dofs[2];
    T width[2];
    T length;

    bool operator<(const PIPE_CHANGE_KEY& p) const
    {
        for(int i=0;i<2;i++)
        {
            if(num_dofs[i]<p.num_dofs[i]) return true;
            if(p.num_dofs[i]<num_dofs[i]) return false;
            if(width[i]<p.width[i]-comp_tol) return true;
            if(p.width[i]<width[i]-comp_tol) return false;
        }
        if(length<p.length-comp_tol) return true;
        return false;
    }
};

template<class T>
struct COMPONENT_CHANGE
{
    typedef VECTOR<T,2> TV;
    
    T target_length;
    std::map<PIPE_CHANGE_KEY<T>,CANONICAL_COMPONENT<T>*> canonical_changes;
    std::map<PIPE_CHANGE_KEY<T>,CANONICAL_BLOCK<T>*> canonical_change_blocks;

    CANONICAL_COMPONENT<T>* Make_Component(const PIPE_CHANGE_KEY<T>& key);
    CANONICAL_BLOCK<T>* Make_Block(const PIPE_CHANGE_KEY<T>& key);
};

}

#endif
