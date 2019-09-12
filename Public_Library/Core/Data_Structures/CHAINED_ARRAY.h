//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CHAINED_ARRAY
//#####################################################################
#ifndef __CHAINED_ARRAY__
#define __CHAINED_ARRAY__

#include <Core/Arrays/ARRAY.h>
namespace PhysBAM{

// T must be:
// (1) an int
// (2) a structure whose first entry is an int, non-virtual.
// The int MUST be non-negative. (the sign bit must never be set.)
//
// Optimized for the case where almost all entries of the array
// contain zero or one element.

template<class T,class ID>
class CHAINED_ARRAY
{
public:
    union DATA
    {
        int i;
        T t;

        DATA(): t() {i=-1;}
        ~DATA() {t.~T();}
    };

    HASHTABLE<ID,DATA> hash;
    ARRAY<ARRAY<T> > extra;

    CHAINED_ARRAY()
    {}

    void Remove_All()
    {hash.Remove_All();extra.Remove_All();}

    template<class... Args>
    void Resize(Args&&... args)
    {Remove_All();}

    template<class U,class CA>
    static ARRAY_VIEW<U> Get_Helper(CA& ca,const ID& i)
    {
        if(DATA* d=ca.hash.Get_Pointer(i)){
            if(d->i>=0) return ARRAY_VIEW<U>(&d->t,1);
            return ca.extra(~d->i);}
        return ARRAY_VIEW<U>();
    }

    ARRAY_VIEW<const T> Get(const ID& i) const
    {return Get_Helper<const T>(*this,i);}
    
    ARRAY_VIEW<T> Get(const ID& i)
    {return Get_Helper<T>(*this,i);}

    void Insert(const ID& i,const T& data)
    {
        if(DATA* d=hash.Get_Pointer(i)){
            if(d->i<0) extra(~d->i).Append(data);
            else{
                int j=extra.Append(ARRAY<T>());
                extra(j).Append(d->t);
                extra(j).Append(data);
                d->i=~j;}}
        else{
            DATA e;
            e.t=data;
            hash.Set(i,e);}
    }

    bool Contains(const ID& i) const
    {return array(i).i!=-1;}

//#####################################################################
};
}
#endif

