//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
// Class CODE_BLOCK_NODE
//#####################################################################
#ifndef __CODE_BLOCK_NODE__
#define __CODE_BLOCK_NODE__

#include <Tools/Symbolics/INSTRUCTION.h>
namespace PhysBAM{

struct CODE_BLOCK;

struct CODE_BLOCK_NODE
{
    CODE_BLOCK_NODE *next,*prev;
    CODE_BLOCK* block;

    INSTRUCTION inst;

    void Set_Use_Def();
    void Unset_Use_Def();

    CODE_BLOCK_NODE()
        :next(0),prev(0),block(0)
    {}
};
}
#endif
