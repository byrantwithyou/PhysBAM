//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COMPONENT_BC__
#define __COMPONENT_BC__

#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <map>
#include "CANONICAL_BLOCK.h"
#include "CANONICAL_COMPONENT.h"
#include "COMPONENT_PIPE.h"

namespace PhysBAM{

template<class T>
struct COMPONENT_BC
{
    typedef VECTOR<T,2> TV;
    typedef TRIPLE<CANONICAL_BLOCK<T>*,INTERVAL<int>,INTERVAL<int> > TRIP;
    
    T target_length;
    std::map<PIPE_KEY<T>,TRIP> canonical_bc_blocks[2];

    TRIP Make_Block(int d,T w,T l,bool is_v);
};

}

#endif
