//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPONENT_PIPE__
#define __COMPONENT_PIPE__

#include <map>
#include "CANONICAL_BLOCK.h"
#include "CANONICAL_COMPONENT.h"

namespace PhysBAM{

extern double comp_tol;

template<class T>
struct PIPE_KEY
{
    int num_dofs;
    T width;
    T length;

    bool operator<(const PIPE_KEY& p) const
    {
        if(num_dofs<p.num_dofs) return true;
        if(p.num_dofs<num_dofs) return false;
        if(width<p.width-comp_tol) return true;
        if(p.width<width-comp_tol) return false;
        if(length<p.length-comp_tol) return true;
        return false;
    }
};

template<class T>
struct COMPONENT_PIPE
{
    typedef VECTOR<T,2> TV;
    
    T target_length;
    std::map<PIPE_KEY<T>,CANONICAL_COMPONENT<T>*> canonical_pipes;
    std::map<PIPE_KEY<T>,CANONICAL_BLOCK<T>*> canonical_pipe_blocks;

    CANONICAL_COMPONENT<T>* Make_Component(const PIPE_KEY<T>& key);
    CANONICAL_BLOCK<T>* Make_Block(const PIPE_KEY<T>& key);
};

}

#endif
