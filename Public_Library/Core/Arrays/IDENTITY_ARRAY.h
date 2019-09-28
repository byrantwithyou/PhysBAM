//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IDENTITY_ARRAY
//#####################################################################
#ifndef __IDENTITY_ARRAY__
#define __IDENTITY_ARRAY__

#include <Core/Arrays/ARRAY_BASE.h>
#include <Core/Arrays/SIMPLE_ITERATOR.h>
#include <cassert>
namespace PhysBAM{

template<class ID> struct IS_ARRAY<IDENTITY_ARRAY<ID> > {static const bool value=true;};
template<class ID> struct IS_ARRAY_VIEW<IDENTITY_ARRAY<ID> > {static const bool value=true;};

template<class ID> // ID=int
class IDENTITY_ARRAY:public ARRAY_BASE<ID,IDENTITY_ARRAY<ID>,ID>
{
public:
    typedef ID ELEMENT;
    typedef ID INDEX;
private:
    ID m;
public:
    typedef SIMPLE_ITERATOR<const IDENTITY_ARRAY> iterator;
    typedef SIMPLE_ITERATOR<const IDENTITY_ARRAY> const_iterator;

    explicit IDENTITY_ARRAY(const ID m)
        :m(m)
    {}

    ID Size() const
    {return m;}

    ID operator()(const ID i) const
    {assert((unsigned)Value(i)<(unsigned)Value(m));return i;}

    struct ITERATOR
    {
        ID i;
        explicit ITERATOR(ID i): i(i) {}
        ID operator*() const {return i;}
        ITERATOR& operator++(){i++;return *this;}
        ITERATOR operator++(int){ITERATOR it(i++);}
        ITERATOR& operator--(){i--;return *this;}
        ITERATOR operator--(int){ITERATOR it(i--);}
        bool operator==(const ITERATOR& it) const {return i==it.i;}
        bool operator!=(const ITERATOR& it) const {return i!=it.i;}
        bool operator<(const ITERATOR& it) const {return i<it.i;}
        bool operator>(const ITERATOR& it) const {return i>it.i;}
        bool operator<=(const ITERATOR& it) const {return i<=it.i;}
        bool operator>=(const ITERATOR& it) const {return i>=it.i;}
    };

    ITERATOR begin() const {return ITERATOR(0);}
    ITERATOR end() const {return ITERATOR(m);}

//#####################################################################
};
}
#endif
