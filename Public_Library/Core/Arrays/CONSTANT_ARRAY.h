//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTANT_ARRAY
//#####################################################################
#ifndef __CONSTANT_ARRAY__
#define __CONSTANT_ARRAY__

#include <Core/Arrays/ARRAY_BASE.h>
namespace PhysBAM{

template<class T,class ID> struct IS_ARRAY<CONSTANT_ARRAY<T,ID> > {static const bool value=true;};
template<class T,class ID> struct IS_ARRAY_VIEW<CONSTANT_ARRAY<T,ID> > {static const bool value=true;};

template<class T,class ID> // ID=int
class CONSTANT_ARRAY:public ARRAY_BASE<T,CONSTANT_ARRAY<T,ID>,ID>
{
public:
    typedef T ELEMENT;
    typedef ID INDEX;
private:
    ID m;
    T constant;
public:
    typedef SIMPLE_ITERATOR<const CONSTANT_ARRAY> iterator;
    typedef SIMPLE_ITERATOR<const CONSTANT_ARRAY> const_iterator;

    CONSTANT_ARRAY(const ID m,const T constant)
        :m(m),constant(constant)
    {}

    ID Size() const
    {return m;}

    const T& operator()(const ID i) const
    {assert((unsigned)Value(i)<(unsigned)Value(m));return constant;}

    SIMPLE_ITERATOR<const CONSTANT_ARRAY> begin() const
    {return SIMPLE_ITERATOR<const CONSTANT_ARRAY>(*this,0);}

    SIMPLE_ITERATOR<const CONSTANT_ARRAY> end() const
    {return SIMPLE_ITERATOR<const CONSTANT_ARRAY>(*this,Size());}

//#####################################################################
};
}
#endif
