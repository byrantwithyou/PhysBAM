//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
// Class CODE_BLOCK
//#####################################################################
#ifndef __CODE_BLOCK__
#define __CODE_BLOCK__

#include <Tools/Data_Structures/HASHTABLE.h>
namespace PhysBAM{

struct CODE_BLOCK_NODE;
struct INSTRUCTION;

struct CODE_BLOCK
{
    CODE_BLOCK *prev[2],*next[2];
    CODE_BLOCK_NODE *head,*tail;
    HASHTABLE<int> in,out;
    int id;

    CODE_BLOCK()
        :head(0),tail(0),id(-1)
    {prev[0]=0;prev[1]=0;next[0]=0;next[1]=0;}

    void Append(CODE_BLOCK_NODE* B);
    void Prepend(CODE_BLOCK_NODE* B);
    void Insert_Before(CODE_BLOCK_NODE* N,CODE_BLOCK_NODE* B); // Insert N before B
    void Insert_After(CODE_BLOCK_NODE* N,CODE_BLOCK_NODE* B); // Insert N after B
    CODE_BLOCK_NODE* Remove_Node(CODE_BLOCK_NODE* N); // Remove from list but do not delete
    void Delete_Node(CODE_BLOCK_NODE* N); // Remove from list and delete
    CODE_BLOCK_NODE* Append(const INSTRUCTION& I);
    CODE_BLOCK_NODE* Prepend(const INSTRUCTION& I);
    CODE_BLOCK_NODE* Insert_Before(const INSTRUCTION& I,CODE_BLOCK_NODE* B); // Insert I before B
    CODE_BLOCK_NODE* Insert_After(const INSTRUCTION& I,CODE_BLOCK_NODE* B); // Insert I after B
    void Merge_Code(CODE_BLOCK* cb);
    void Check_Links() const;
};
}
#endif
