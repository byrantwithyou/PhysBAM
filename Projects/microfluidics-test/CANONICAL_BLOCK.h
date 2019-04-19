//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CANONICAL_BLOCK__
#define __CANONICAL_BLOCK__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Vectors/VECTOR.h>
#include "COMMON.h"

namespace PhysBAM{

// v = range of vertex indices in cross section
// e = range of edge indices in cross section
// if own_first, first half of v and e is owned by this block
// if v or e has odd size; middle is master
struct CROSS_SECTION
{
    INTERVAL<int> v,e;
    bool own_first;
};

template<class T>
struct CANONICAL_BLOCK
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,3> IV3;
    typedef VECTOR<int,2> IV;

    ARRAY<CROSS_SECTION,CON_ID> cross_sections;
    ARRAY<TV> X;
    ARRAY<IV3> E;
    ARRAY<IV> S;
    ARRAY<int> bc_v,bc_e;
    ARRAY<int> ticks; // edge index e -> t. tick is on the side of S(e)(t)
};
}
#endif
