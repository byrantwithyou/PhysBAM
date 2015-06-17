//#####################################################################
// Copyright 2003-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HASHTABLE_ITERATOR
//##################################################################### 
//
// Use this to iterate through all values stored in the hashtable.
// Loop through using something like:
//     for(HASHTABLE_ITERATOR<TK,T> iterator(hash);iterator.Valid();iterator.Next()){
//         TK key=iterator.Key();T data=iterator.Data();}
//
//#####################################################################
#ifndef __HASHTABLE_ITERATOR__
#define __HASHTABLE_ITERATOR__

#include <Tools/Data_Structures/HASHTABLE.h>
namespace PhysBAM{

template<class TK,class T> // T = void
class HASHTABLE_ITERATOR
{
    struct UNUSABLE{};
    typedef typename remove_const<T>::type RAW_T;
    typedef typename IF<is_same<RAW_T,void>::value,UNUSABLE,T>::TYPE T_UNLESS_VOID;
    typedef typename IF<is_same<RAW_T,void>::value || is_const<T>::value,const HASHTABLE<TK,RAW_T>,HASHTABLE<TK,T> >::TYPE T_HASHTABLE;
public:
    T_HASHTABLE& hashtable;
private:
    int current_h;
public:
    enum INVALID_ENUM {INVALID_ITERATOR};

    HASHTABLE_ITERATOR(T_HASHTABLE& hashtable)
        :hashtable(hashtable)
    {
        Reset();
    }

    HASHTABLE_ITERATOR(T_HASHTABLE& hashtable,INVALID_ENUM)
        :hashtable(hashtable),current_h(-1)
    {}

    void Reset()
    {for(current_h=0;current_h<hashtable.table.m;current_h++) if(hashtable.table(current_h).state==ENTRY_ACTIVE) return;current_h=-1;}

    const TK& Key() const
    {return hashtable.table(current_h).key;}

    T_UNLESS_VOID& Data()
    {return hashtable.table(current_h).data;}

    const T_UNLESS_VOID& Data() const
    {return hashtable.table(current_h).data;}

    bool Valid() const
    {return current_h!=-1;}

    void Next()
    {assert(Valid());current_h++;
    for(;current_h<hashtable.table.m;current_h++) if(hashtable.table(current_h).state==ENTRY_ACTIVE) return;current_h=-1;}

//#####################################################################
};
}
#endif
