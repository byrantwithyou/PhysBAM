//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYREGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Symbolics/CODE_BLOCK.h>
#include <Tools/Symbolics/CODE_BLOCK_NODE.h>
using namespace PhysBAM;
//#####################################################################
// Function Append
//#####################################################################
void CODE_BLOCK::
Append(CODE_BLOCK_NODE* N)
{
    N->block=this;
    if(!tail){
        N->next=N->prev=0;
        tail=head=N;}
    else{
        tail->next=N;
        N->prev=tail;
        N->next=0;
        tail=N;}
}
//#####################################################################
// Function Prepend
//#####################################################################
void CODE_BLOCK::
Prepend(CODE_BLOCK_NODE* N)
{
    N->block=this;
    if(!head){
        N->next=N->prev=0;
        tail=head=N;}
    else{
        head->prev=N;
        N->next=head;
        N->prev=0;
        head=N;}
}
//#####################################################################
// Function Insert_Before
//#####################################################################
// Insert N before B
void CODE_BLOCK::
Insert_Before(CODE_BLOCK_NODE* N,CODE_BLOCK_NODE* B)
{
    N->block=this;
    N->prev=B->prev;
    N->next=B;
    B->prev=N;
    if(N->prev) N->prev->next=N;
    else head=N;
}
//#####################################################################
// Function Insert_After
//#####################################################################
// Insert N after B
void CODE_BLOCK::
Insert_After(CODE_BLOCK_NODE* N,CODE_BLOCK_NODE* B)
{
    N->block=this;
    N->next=B->next;
    N->prev=B;
    B->next=N;
    if(N->next) N->next->prev=N;
    else tail=N;
}
//#####################################################################
// Function Remove_Node
//#####################################################################
// Remove from list but do not delete
CODE_BLOCK_NODE* CODE_BLOCK::
Remove_Node(CODE_BLOCK_NODE* N)
{
    if(N->next) N->next->prev=N->prev;
    else tail=N->prev;
    if(N->prev) N->prev->next=N->next;
    else head=N->next;
    N->block=0;
    N->next=0;
    N->prev=0;
    return N;
}
//#####################################################################
// Function Remove_Node
//#####################################################################
void CODE_BLOCK::
Delete_Node(CODE_BLOCK_NODE* N)
{
    if(N->next) N->next->prev=N->prev;
    else tail=N->prev;
    if(N->prev) N->prev->next=N->next;
    else head=N->next;
    delete N;
}
//#####################################################################
// Function Append
//#####################################################################
CODE_BLOCK_NODE* CODE_BLOCK::
Append(const INSTRUCTION& I)
{
    CODE_BLOCK_NODE* N=new CODE_BLOCK_NODE;
    N->inst=I;
    Append(N);
    return N;
}
//#####################################################################
// Function Prepend
//#####################################################################
CODE_BLOCK_NODE* CODE_BLOCK::
Prepend(const INSTRUCTION& I)
{
    CODE_BLOCK_NODE* N=new CODE_BLOCK_NODE;
    N->inst=I;
    Prepend(N);
    return N;
}
//#####################################################################
// Function Insert_Before
//#####################################################################
// Insert I before B
CODE_BLOCK_NODE* CODE_BLOCK::
Insert_Before(const INSTRUCTION& I,CODE_BLOCK_NODE* B)
{
    CODE_BLOCK_NODE* N=new CODE_BLOCK_NODE;
    N->inst=I;
    Insert_Before(N,B);
    return N;
}
//#####################################################################
// Function Insert_After
//#####################################################################
// Insert I after B
CODE_BLOCK_NODE* CODE_BLOCK::
Insert_After(const INSTRUCTION& I,CODE_BLOCK_NODE* B)
{
    CODE_BLOCK_NODE* N=new CODE_BLOCK_NODE;
    N->inst=I;
    Insert_After(N,B);
    return N;
}
//#####################################################################
// Function Merge_Code
//#####################################################################
void CODE_BLOCK::
Merge_Code(CODE_BLOCK* cb)
{
    if(!cb->head) return;
    if(!head){
        head=cb->head;
        tail=cb->tail;}
    else{
        tail->next=cb->head;
        cb->head->prev=tail;
        tail=cb->tail;}
    for(CODE_BLOCK_NODE* N=cb->head;N;N=N->next) N->block=this;
    cb->head=0;
    cb->tail=0;
}
//#####################################################################
// Function Check_Links
//#####################################################################
void CODE_BLOCK::
Check_Links() const
{
    if(!head){PHYSBAM_ASSERT(!tail);return;}
    PHYSBAM_ASSERT(tail);
    PHYSBAM_ASSERT(tail->next==0);
    PHYSBAM_ASSERT(head->prev==0);

    for(CODE_BLOCK_NODE* N=head;N;N=N->next)
    {
        PHYSBAM_ASSERT(N==(N->next?N->next->prev:tail));
        PHYSBAM_ASSERT(N==(N->prev?N->prev->next:head));
        PHYSBAM_ASSERT(N->block==this);
    }
}
