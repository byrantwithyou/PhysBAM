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

template<class T,class ID=int> // ID=int
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

    ARRAY<DATA,ID> array;
    ARRAY<ARRAY<T> > extra;

    CHAINED_ARRAY()
    {}

    void Remove_All()
    {for(auto& it:array) it.i=-1;}

    template<class... Args>
    void Resize(Args&&... args)
    {array.Resize(args...);Remove_All();}

    template<class U,class CA>
    static ARRAY_VIEW<U> Get_Helper(CA& ca,const ID& i)
    {
        auto& d=ca.array(i);
        if(d.i==-1) return ARRAY_VIEW<U>();
        if(d.i>=0) return ARRAY_VIEW<U>(1,&d.t);
        return ca.extra(-2-d.i);
    }
    
    ARRAY_VIEW<const T> Get(const ID& i) const
    {return Get_Helper<const T>(*this,i);}
    
    ARRAY_VIEW<T> Get(const ID& i)
    {return Get_Helper<T>(*this,i);}

    void Insert(const ID& i,const T& data)
    {
        DATA& d=array(i);
        if(d.i==-1){d.t=data;assert(d.i>=0);}
        else if(d.i<0) extra(-2-d.i).Append(data);
        else{
            int j=extra.Append(ARRAY<T>());
            extra(j).Append(d.t);
            extra(j).Append(data);
            d.i=-2-j;}
    }

    bool Contains(const ID& i) const
    {return array(i).i!=-1;}
    
//#####################################################################
};
}
#endif

