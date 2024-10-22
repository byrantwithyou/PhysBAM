//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace ARRAYS_COMPUTATIONS
//#####################################################################
#ifndef __ARRAY_HEAPIFY__
#define __ARRAY_HEAPIFY__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/exchange.h>
namespace PhysBAM{

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T2,class ID>
    void Heapify(ARRAY<T,ID>& a,ARRAY<T2,ID>& aux,ID index,const ID heap_size) // largest on top, only sorts down from index (not up!)
    {for(;;){ID left(2*Value(index)+1),right(2*Value(index)+2),index_of_largest=index;
        if(left<heap_size && a(left)>a(index_of_largest)) index_of_largest=left;
        if(right<heap_size && a(right)>a(index_of_largest)) index_of_largest=right;
        if(index_of_largest!=index){exchange(a(index),a(index_of_largest));exchange(aux(index),aux(index_of_largest));index=index_of_largest;}else return;}}

    template<class T,class T2,class ID>
    void Heapify(ARRAY<T,ID>& a,ARRAY<T2,ID>& aux) // largest on top
    {for(ID i=a.m/2-1;i>=ID(0);i--) Heapify(a,aux,i,a.m);}

    template<class T,class T2,class ID>
    void Heapify(ARRAY<T,ID>& a,ARRAY<T2,ID>& aux,const ID max_index) // largest on top, only does from 1 to max_index
    {for(ID i(Value((max_index-1)/2));i>=ID(0);i--) Heapify(a,aux,i,max_index);}

    template<class T_ARRAY>
    void Heapify(T_ARRAY& a,typename T_ARRAY::INDEX index,const typename T_ARRAY::INDEX heap_size) // largest on top, only sorts down from index (not up!)
    {typedef typename T_ARRAY::INDEX ID;
    for(;;){ID left(2*Value(index)+1),right(2*Value(index)+2),index_of_largest=index;
        if(left<heap_size && a(left)>a(index_of_largest)) index_of_largest=left;
        if(right<heap_size && a(right)>a(index_of_largest)) index_of_largest=right;
        if(index_of_largest!=index){exchange(a(index),a(index_of_largest));index=index_of_largest;}else return;}}

    template<class T_ARRAY>
    void Heapify(T_ARRAY& a) // largest on top
    {typedef typename T_ARRAY::INDEX ID;for(ID i=a.Size()/2-1;i>=ID(0);i--) Heapify(a,i,a.Size());}

    template<class T,class ID>
    void Heapify(ARRAY<T,ID>& a,const ID max_index) // largest on top, only does from 0 to max_index-1
    {for(ID i=max_index/2-1;i>=ID(0);i--) Heapify(a,i,max_index);}

//#####################################################################
}
}
#endif
