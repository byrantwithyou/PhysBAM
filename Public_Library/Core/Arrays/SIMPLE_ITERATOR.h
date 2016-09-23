//#####################################################################
// Copyright 2004-2009, Kevin Der, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_ITERATOR
//#####################################################################
#ifndef __SIMPLE_ITERATOR__
#define __SIMPLE_ITERATOR__
namespace PhysBAM{

template<class T_ARRAY>
struct SIMPLE_ITERATOR
{
    T_ARRAY& array;
    typename T_ARRAY::INDEX index;

    SIMPLE_ITERATOR(T_ARRAY& array,typename T_ARRAY::INDEX index)
        :array(array),index(index)
    {}

    typename T_ARRAY::ELEMENT& operator*()
    {return array(index);}

    const typename T_ARRAY::ELEMENT& operator*() const
    {return array(index);}

    typename T_ARRAY::ELEMENT* operator->()
    {return &array(index);}

    const typename T_ARRAY::ELEMENT* operator->() const
    {return &array(index);}

    // stl
    SIMPLE_ITERATOR& operator++()
    {index++;return *this;}

    SIMPLE_ITERATOR operator++(int)
    {SIMPLE_ITERATOR it(*this);index++;return it;}

    // stl
    SIMPLE_ITERATOR& operator--()
    {index--;return *this;}

    SIMPLE_ITERATOR operator--(int)
    {SIMPLE_ITERATOR it(*this);index--;return it;}

    template<class T_ARRAY1> bool operator==(const SIMPLE_ITERATOR<T_ARRAY1>& it) const
    {return index==it.index;}

    template<class T_ARRAY1> bool operator!=(const SIMPLE_ITERATOR<T_ARRAY1>& it) const
    {return index!=it.index;}

    template<class T_ARRAY1> bool operator<(const SIMPLE_ITERATOR<T_ARRAY1>& it) const
    {return index<it.index;}

    template<class T_ARRAY1> bool operator>(const SIMPLE_ITERATOR<T_ARRAY1>& it) const
    {return index>it.index;}

    template<class T_ARRAY1> bool operator<=(const SIMPLE_ITERATOR<T_ARRAY1>& it) const
    {return index<=it.index;}

    template<class T_ARRAY1> bool operator>=(const SIMPLE_ITERATOR<T_ARRAY1>& it) const
    {return index>=it.index;}

};
}
#endif
