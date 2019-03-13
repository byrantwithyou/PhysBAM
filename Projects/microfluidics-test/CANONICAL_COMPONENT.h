//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CANONICAL_COMPONENT__
#define __CANONICAL_COMPONENT__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/ELEMENT_ID.h>
#include <Core/Data_Structures/PAIR.h>
#include <map>
#include "COMMON.h"
#include "XFORM.h"
namespace PhysBAM{

PHYSBAM_DECLARE_ELEMENT_ID(CC_BLOCK_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(CC_IRREG_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);

template<class T> struct CANONICAL_BLOCK;

// first is master, second is slave
struct CC_BLOCK_CONNECTION
{
    CC_BLOCK_ID id;
    CON_ID con_id;
    CC_IRREG_ID irreg_id;
    bool master,is_regular;

    CC_BLOCK_CONNECTION(CC_BLOCK_ID b=CC_BLOCK_ID(-7),CON_ID c=CON_ID(-7))
        :id(b),con_id(c),irreg_id(-7),master(false),is_regular(true)
    {}

    CC_BLOCK_CONNECTION(CC_IRREG_ID i)
        :id(-7),con_id(-7),irreg_id(i),master(false),is_regular(false)
    {}

    void Set_Irreg(CC_IRREG_ID i)
    {
        con_id=CON_ID(-7);
        irreg_id=i;
        is_regular=false;
    }
};

template<class T>
struct CC_BLOCK
{
    typedef VECTOR<T,2> TV;
    CANONICAL_BLOCK<T>* block;
    XFORM<TV> xform;
    ARRAY<CC_BLOCK_CONNECTION,CON_ID> connections;
    ARRAY<CC_IRREG_ID> edge_on; // for edge-on (index in irregular_connections)
    int flags=0; // 1=separator, 2=separator-eligible
};

// regular is master
struct CC_IRREGULAR_CONNECTION
{
    CC_BLOCK_ID regular=CC_BLOCK_ID(-7);
    CON_ID con_id;
    // one for each dof on cross section, starting from owned side of cross section
    ARRAY<PAIR<CC_BLOCK_ID,int> > edge_on_v,edge_on_e;
    CC_IRREG_ID ref_ic=CC_IRREG_ID(-7);
};

// neighbor block i is given index ~i and con_id=-1.
template<class T>
struct CANONICAL_COMPONENT
{
    ARRAY<CC_BLOCK<T>,CC_BLOCK_ID> blocks;
    ARRAY<CC_IRREGULAR_CONNECTION,CC_IRREG_ID> irregular_connections;
};
}
#endif
